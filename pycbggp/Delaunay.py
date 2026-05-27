
from Geometry import Point as P

arb = P(10**18, 10**18)

class Quad:
  def __init__(self, s = None):
    self.rot = s 
    self.o = None
    self.p = arb
    self.mark = False 

  def r(self) -> 'Quad':
    return self.rot.rot
  def set_r(self, q):
    self.rot.rot = q
  def F(self) -> P:
    return self.r().p
  def set_F(self, p : P):
    self.r().p = p
  def prev(self) -> 'Quad':
    return self.rot.o.rot
  def set_prev(self, q):
    self.rot.o.rot = q
  def next(self) -> 'Quad':
    return self.r().prev()
  def set_next(self, q):
    self.r().set_prev(q)
  
class Delaunay:
  def __init__(self, points = None):
    # self.points = points
    self.H = None
  
  # is p in the circumcircle of a, b, c?
  def circ(self, p : P, a : P, b : P, c : P):
    p2 = p.dist2()
    A = a.dist2() - p2
    B = b.dist2() - p2 
    C = c.dist2() - p2
    return p.cross(a, b) * C + p.cross(b, c) * A + p.cross(c, a) * B > 0

  def make_edge(self, a : P, b : P):
    r = self.H if self.H  else Quad(Quad(Quad(Quad())))
    self.H = r.o
    # r.r().r() = r
    r.r().set_r(r)
    for i in range(4):
      r = r.rot
      r.p = arb
      r.o = r if i & 1 == 1 else r.r()

    r.p = a
    # r.F() = b
    r.set_F(b)
    return r
  
  def splice(self, a : Quad, b : Quad):
    a.o.rot.o, b.o.rot.o = b.o.rot.o, a.o.rot.o
    a.o, b.o = b.o, a.o
  
  def connect(self, a : Quad, b : Quad):
    q = self.make_edge(a.F(), b.p)
    self.splice(q, a.next())
    self.splice(q.r(), b)
    return q
  
  def rec(self, s : list[P]):
    if len(s) < 4:
      a = self.make_edge(s[0], s[1])
      if (len(s) == 2):
        return (a, a.r())
      b = self.make_edge(s[1], s[-1])
      self.splice(a.r(), b)
      side = s[0].cross(s[1], s[2])
      c = self.connect(b, a) if side else None
      return ((c.r() if side < 0 else a), (c if side < 0 else b.r()))
    
    half = (len(s) + 1) >> 1
    ra, A = self.rec(s[:half])
    B, rb = self.rec(s[half:])

    while True:
      if (B.p.cross(A.F(), A.p) < 0):
        A = A.next()
      elif (A.p.cross(B.F(), B.p) > 0):
        B = B.r().o
      else:
        break

    base = self.connect(B.r(), A)
    if (A.p == ra.p):
      ra = base.r()
    if (B.p == rb.p):
      rb = base
    
    while True:
      LC = base.r().o
      if LC.F().cross(base.F(), base.p) > 0:
        while True:
          if (not self.circ(LC.o.F(), base.F(), base.p, LC.F())):
            break
          t = LC.o
          self.splice(LC, LC.prev())
          self.splice(LC.r(), LC.r().prev())
          LC.o = self.H
          self.H = LC
          LC = t

      RC = base.prev()
      if RC.F().cross(base.F(), base.p) > 0:
        while True:
          if (not self.circ(RC.prev().F(), base.F(), base.p, RC.F())):
            break
          t = RC.prev()
          self.splice(RC, RC.prev())
          self.splice(RC.r(), RC.r().prev())
          RC.o = self.H
          self.H = RC
          RC = t
      
      if (LC.F().cross(base.F(), base.p) <= 0 and RC.F().cross(base.F(), base.p) <= 0):
        break
      if (LC.F().cross(base.F(), base.p) <= 0 or (RC.F().cross(base.F(), base.p) > 0 and self.circ(RC.F(), RC.p, LC.F(), LC.p))):
        base = self.connect(RC, base.r())
      else:
        base = self.connect(base.r(), LC.r())

    return (ra, rb)
  
  def triangulate(self, pts : list[P]):
    pts = sorted(pts)
    self.H = None
    assert all(pts[i] != pts[i - 1] for i in range(1, len(pts)))
    if (len(pts) < 2):
      pts.clear()
      return pts
    (e, _) = self.rec(pts)
    q = [e]

    qi = 0
    while e.o.F().cross(e.F(), e.p) < 0:
      e = e.o
    
    c = e
    while True:
      c.mark = True
      pts.append(c.p)
      q.append(c.r())
      c = c.next()
      if (c == e):
        break

    pts.clear()
    while qi < len(q):
      e = q[qi]
      if (not e.mark):
        c = e
        while True:
          c.mark = True
          pts.append(c.p)
          q.append(c.r())
          c = c.next()
          if (c == e):
            break
      qi += 1

    return pts






