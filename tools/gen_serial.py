#!/usr/bin/env python3
"""Generate pic_mpi_serial from pic_mpi_f08's signatures.

The serial backend has to present the SAME interface as the MPI one -- 140
procedures behind a dozen generics -- and only the bodies differ. Retyping
those signatures by hand would be 140 chances to introduce a mismatch that
only shows up as an unresolved generic in somebody's build.

So the signatures are taken from the real module and only the bodies are
written here, by family.
"""
import re
import sys

SRC = sys.argv[1]
OUT = sys.argv[2]

lines = open(SRC).read().splitlines()
#
# The MODULE-level `contains`, not the one inside a derived type.
#
# `type :: comm_t` has its own `contains` for the type-bound procedures, and
# it comes first in the file -- so taking the first match copies half a type
# definition and leaves `comm_t` undefined.
#
depth = 0
ci = None
for i, l in enumerate(lines):
    s_ = l.strip().lower()
    if re.match(r'type\s*(,.*)?::\s*\w+', s_) and 'end type' not in s_:
        depth += 1
    elif s_.startswith('end type'):
        depth -= 1
    elif s_ == 'contains' and depth == 0:
        ci = i
        break
if ci is None:
    raise SystemExit('no module-level contains found')

PROC_RE = re.compile(r'\s*(pure\s+|elemental\s+)?(subroutine|function)\s+(\w+)', re.I)
DECL_RE = re.compile(
    r'^\s*($|!|&|use\b|implicit\b|import\b|'
    r'(type|class|integer|real|logical|character|complex|procedure)\s*[\(,:]|'
    r'(type|class)\s*\(.*\)\s*(,|::))', re.I)

procs = []
i = ci + 1
while i < len(lines):
    m = PROC_RE.match(lines[i])
    if m:
        kind, name = m.group(2).lower(), m.group(3)
        j, depth = i, 0
        while j < len(lines):
            s = lines[j].strip().lower()
            if PROC_RE.match(lines[j]) and not s.startswith('end'):
                depth += 1
            if re.match(r'end\s*(subroutine|function)\b', s):
                depth -= 1
                if depth == 0:
                    break
            j += 1
        # header runs while lines continue with `&`
        h = i
        while lines[h].rstrip().endswith('&'):
            h += 1
        header = lines[i:h + 1]
        # declarations: from after the header while lines look declarative
        d = h + 1
        decls = []
        while d < j:
            if DECL_RE.match(lines[d]):
                # A procedure-local `use mpi_f08, only: ...` has nothing to
                # resolve against in this module; the body that needed it is
                # being replaced anyway.
                if not re.match(r'\s*use\s+mpi', lines[d], re.I):
                    decls.append(lines[d])
                # carry continuations
                while lines[d].rstrip().endswith('&'):
                    d += 1
                    decls.append(lines[d])
                d += 1
            else:
                break
        result = None
        rm = re.search(r'result\s*\(\s*(\w+)\s*\)', ' '.join(header), re.I)
        if rm:
            result = rm.group(1)
        elif kind == 'function':
            result = name
        procs.append(dict(name=name, kind=kind, header=header, decls=decls,
                          result=result, end=lines[j]))
        i = j + 1
        continue
    i += 1

def dummies(header):
    txt = ' '.join(h.replace('&', ' ') for h in header)
    m = re.search(r'\((.*?)\)', txt)
    return [a.strip() for a in m.group(1).split(',')] if m else []

def abort(n):
    """A stop that names the call and says why it cannot be reached."""
    return [
        '      ! A single rank has no peer. Reaching this means the caller',
        '      ! did not take its `size() == 1` path -- which is a bug worth',
        '      ! stopping for, not one to hide by returning quietly.',
        f'      error stop "pic_mpi_serial: {n} needs a peer rank"',
    ]

def body(p):
    n, out = p['name'], []
    dl = dummies(p['header'])
    def setres(v):
        return [f'      {p["result"]} = {v}'] if p['result'] else []

    # --- point to point and one-sided: unreachable at one rank ------------
    if re.match(r'comm_(send|isend|recv|irecv)_', n) or n.startswith('win_') \
            or n in ('comm_probe',):
        if n in ('win_is_null',):
            return setres('.true.')
        if n == 'win_get_handle':
            return setres('MPI_WIN_NULL')
        return abort(n)

    # --- collectives: on one rank the answer is already in place ----------
    if re.match(r'comm_bcast_', n):
        return ['      ! One rank: the data is already where it needs to be.']
    if re.match(r'comm_(allgather|allreduce)_', n):
        # copy send buffer into recv buffer if both are present
        cand = [a for a in dl if a not in ('comm', 'op', 'root', 'count')]
        if len(cand) >= 2:
            s_, r_ = cand[0], cand[1]
            return [f'      ! One rank: the gathered/reduced result is this',
                    f'      ! rank\'s own contribution.',
                    f'      {r_} = {s_}']
        return ['      ! One rank: nothing to combine.']

    # --- probes and completion -------------------------------------------
    if n.startswith('comm_iprobe'):
        # No peer can ever have sent us anything. Find the logical out-argument
        # from the declarations rather than guessing what it is called.
        for d in p['decls']:
            m = re.match(r'\s*logical[^:]*intent\s*\(\s*out\s*\)[^:]*::\s*(\w+)',
                         d, re.I)
            if m:
                return ['      ! No peer exists, so nothing is ever pending.',
                        '      ! Reporting false keeps a polling loop terminating',
                        '      ! rather than spinning forever.',
                        f'      {m.group(1)} = .false.']
        return setres('.false.')
    if n.startswith('request_wait') or n == 'wait_all' or 'waitall' in n:
        return ['      ! Nothing was ever started, so nothing is outstanding.']
    if n.startswith('request_test'):
        for d in p['decls']:
            m = re.match(r'\s*logical[^:]*intent\s*\(\s*out\s*\)[^:]*::\s*(\w+)',
                         d, re.I)
            if m:
                return ['      ! Nothing was started, so it is complete.',
                        f'      {m.group(1)} = .true.']
        return setres('.true.')

    # --- communicator queries --------------------------------------------
    simple = {
        'comm_rank': '0', 'm_size_func': '1', 'comm_leader': '.true.',
        'comm_is_null': '.not. this%is_valid',
    }
    if n in simple:
        return setres(simple[n])
    if n in ('comm_barrier', 'comm_finalize', 'pic_mpi_finalize'):
        return ['      ! Nothing to synchronise or release on one rank.']
    if n == 'pic_mpi_init':
        out = ['      ! No MPI to start. A serial build is always as',
               '      ! thread-safe as it is ever going to be.']
        for d in p['decls']:
            m = re.match(r'\s*integer[^:]*intent\s*\(\s*out\s*\)[^:]*::\s*(\w+)',
                         d, re.I)
            if m:
                opt = 'optional' in d.lower()
                v = m.group(1)
                out.append(f'      if (present({v})) {v} = MPI_THREAD_MULTIPLE'
                           if opt else f'      {v} = MPI_THREAD_MULTIPLE')
        return out
    if n == 'pic_mpi_query_thread_level':
        return setres('MPI_THREAD_MULTIPLE')
    if n == 'abort_comm':
        return ['      ! No ranks to bring down with us.',
                '      error stop "pic_mpi_serial: abort_comm"']
    if n == 'get_processor_name':
        for a in dl:
            if 'name' in a.lower():
                return [f'      {a} = "serial"']
    if n in ('create_world_comm', 'create_null_comm', 'create_comm_from_mpi'):
        return ['      ! The only communicator there is.',
                f'      {p["result"] or "comm"}%m_rank = 0',
                f'      {p["result"] or "comm"}%m_size = 1',
                f'      {p["result"] or "comm"}%is_valid = '
                + ('.false.' if 'null' in n else '.true.')]
    if n.startswith('comm_split') or n.startswith('comm_discard') \
            or n == 'comm_duplicate':
        return ['      ! Every split of one rank is that rank.',
                f'      {p["result"] or "new_comm"} = this']
    if n == 'comm_get':
        return setres('this%m_comm')

    # --- anything unclassified: say so loudly rather than guess ----------
    return abort(n)

o = []
o.append('''!> Single-rank stand-in for `pic_mpi_f08`, for builds without MPI.
!!
!! WHY THIS EXISTS
!! ---------------
!! Not every consumer of PIC wants an MPI dependency. A CI runner without one,
!! a laptop build, or a compiler whose vendor MPI is awkward to get hold of --
!! in each case the alternative is no build at all, because `pic_mpi_lib` had
!! only the two MPI backends to choose from.
!!
!! This is the third: the same public interface, with one rank and no library
!! underneath. `pic_mpi_lib` selects it with `USE_SERIAL`, exactly as it
!! selects the legacy bindings with `USE_LEGACY`, so a consumer's source does
!! not change.
!!
!! WHAT THE BODIES DO
!! ------------------
!!   * `rank` is 0, `size` is 1, and every rank is the leader.
!!   * A broadcast is a no-op: the data is already on the only rank there is.
!!   * A gather or a reduction copies the contribution into the result.
!!   * `iprobe` reports nothing pending, because no peer exists to have sent
!!     anything. That keeps polling loops terminating rather than spinning.
!!   * Waiting completes immediately: nothing was ever started.
!!   * POINT-TO-POINT AND ONE-SIDED CALLS STOP THE PROGRAM. A single rank has
!!     no peer, so reaching a send or a get means the caller failed to take
!!     its `size() == 1` path. Returning quietly would turn that into wrong
!!     results much later; stopping names it where it happened.
!!
!! GENERATED, and deliberately. The signatures are taken from
!! `mpi_f08/pic_mpi.f90` so that the two cannot disagree -- 140 procedures
!! behind a dozen generics is 140 chances to typo an interface that would only
!! surface as an unresolved generic in somebody else's build. Only the bodies
!! are written by hand, by family. Regenerate with `tools/gen_serial.py`.
!!
module pic_mpi_serial
   use pic_types, only: int32, int64, sp, dp
   implicit none
   private
''')
# public list, copied from the real module
pub = [l for l in lines[:ci] if re.match(r'\s*public\s*::', l)]
o.extend(pub)
o.append('''
   !
   ! Stand-ins for the handles and constants `mpi_f08` would supply. They carry
   ! no meaning here beyond letting the same declarations compile: a consumer
   ! that stores an `MPI_Status` or compares against `MPI_ANY_SOURCE` keeps
   ! working, it just never receives anything.
   !
   type :: MPI_Comm
      integer(int32) :: v = 0
   end type MPI_Comm

   type :: MPI_Request
      integer(int32) :: v = 0
   end type MPI_Request

   type :: MPI_Status
      integer(int32) :: MPI_SOURCE = -1
      integer(int32) :: MPI_TAG = -1
      integer(int32) :: MPI_ERROR = 0
   end type MPI_Status

   type :: MPI_Op
      integer(int32) :: v = 0
   end type MPI_Op

   type :: MPI_Win
      integer(int32) :: v = 0
   end type MPI_Win

   type(MPI_Comm), parameter :: MPI_COMM_NULL = MPI_Comm(0)
   type(MPI_Op), parameter :: MPI_SUM = MPI_Op(1)
   type(MPI_Op), parameter :: MPI_MIN = MPI_Op(2)
   type(MPI_Op), parameter :: MPI_MAX = MPI_Op(3)

   integer(int32), parameter :: MPI_ANY_SOURCE = -1
   integer(int32), parameter :: MPI_ANY_TAG = -1
   integer(int32), parameter :: MPI_MAX_PROCESSOR_NAME = 256
   integer(int32), parameter :: MPI_THREAD_SINGLE = 0
   integer(int32), parameter :: MPI_THREAD_FUNNELED = 1
   integer(int32), parameter :: MPI_THREAD_SERIALIZED = 2
   integer(int32), parameter :: MPI_THREAD_MULTIPLE = 3
   integer, parameter :: MPI_ADDRESS_KIND = int64

   type(MPI_Request), parameter :: MPI_REQUEST_NULL = MPI_Request(0)
   type(MPI_Win), parameter :: MPI_WIN_NULL = MPI_Win(0)
''')
# types and interfaces, verbatim from the real module
tstart = next(i for i, l in enumerate(lines) if re.match(r'\s*type\s*::\s*request_t', l))
o.extend(l for l in lines[tstart:ci]
         if not re.match(r'\s*use\s+mpi', l, re.I))
o.append('')
o.append('contains')
empty = []
for p in procs:
    b = body(p)
    if not [x for x in b if x.strip() and not x.strip().startswith('!')]:
        if not re.match(r'comm_bcast_|comm_barrier|comm_finalize|pic_mpi_finalize|request_wait|wait_all|waitall', p['name']):
            empty.append(p['name'])
if empty:
    print('EMPTY BODIES (classify these):', ', '.join(empty))

for p in procs:
    o.append('')
    o.extend(p['header'])
    o.extend(p['decls'])
    o.extend(body(p))
    o.append(p['end'])
o.append('')
o.append('end module pic_mpi_serial')

open(OUT, 'w').write('\n'.join(o) + '\n')
print(f'wrote {OUT}: {len(procs)} procedures, {len(o)} lines')
