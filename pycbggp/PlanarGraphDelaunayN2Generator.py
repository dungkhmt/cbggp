from CBSGG import Graph
import random
import math
from typing import List, Tuple, Optional

Point = Tuple[float, float]
Edge = Tuple[int, int]
Tri = Tuple[int, int, int]

EPS = 1e-12


def circumcircle(ax: float, ay: float, bx: float, by: float, cx: float, cy: float) -> Tuple[float, float, float]:
  d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
  if abs(d) < 1e-18:
    return 0.0, 0.0, float("inf")
  a2 = ax * ax + ay * ay
  b2 = bx * bx + by * by
  c2 = cx * cx + cy * cy
  ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d
  uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d
  r2 = (ux - ax) ** 2 + (uy - ay) ** 2
  return ux, uy, r2


def tri_edges(a: int, b: int, c: int) -> List[Edge]:
  e = [(a, b), (b, c), (c, a)]
  return [tuple(sorted(x)) for x in e]


def delaunay_bowyer_watson(points: List[Point]) -> List[Tri]:
  n = len(points)
  if n < 3:
    return []

  xs = [p[0] for p in points]
  ys = [p[1] for p in points]
  minx, maxx = min(xs), max(xs)
  miny, maxy = min(ys), max(ys)
  dx, dy = maxx - minx, maxy - miny
  D = max(dx, dy)
  cx, cy = (minx + maxx)/2.0, (miny + maxy)/2.0

  st = [(cx, cy + 3*D + 10.0),
        (cx - 3*D - 10.0, cy - 3*D - 10.0),
        (cx + 3*D + 10.0, cy - 3*D - 10.0)]

  P = points + st
  s0, s1, s2 = n, n+1, n+2

  tris: List[Tri] = [(s0, s1, s2)]

  for pi in range(n):
    px, py = P[pi]
    bad: List[Tri] = []
    for (a, b, c) in tris:
      ux, uy, r2 = circumcircle(P[a][0], P[a][1],
                                P[b][0], P[b][1],
                                P[c][0], P[c][1])
      if r2 == float("inf"):
        continue
      if (px - ux)**2 + (py - uy)**2 <= r2 + 1e-12:
        bad.append((a, b, c))

    edge_count: dict[Edge, int] = {}

    def add_edge(u: int, v: int) -> None:
      if u > v:
        u, v = v, u
      edge_count[(u, v)] = edge_count.get((u, v), 0) + 1

    for (a, b, c) in bad:
      add_edge(a, b)
      add_edge(b, c)
      add_edge(c, a)

    bad_set = set(bad)
    tris = [t for t in tris if t not in bad_set]

    boundary = [e for e, cnt in edge_count.items() if cnt == 1]

    for (u, v) in boundary:
      tris.append((u, v, pi))

  tris = [t for t in tris if t[0] < n and t[1] < n and t[2] < n]
  return tris


def edges_from_tris(tris: List[Tri]) -> List[Edge]:
  S = set()
  for a, b, c in tris:
    for e in tri_edges(a, b, c):
      S.add(e)
  return sorted(S)


class PlanarGraphTriangularGenerator:
  # generate planar graph n nodes, m edges, max degree d
  # triangular method
  def Generate(self, n, m, d, seed=998244353):
    if m > 3 * n - 6:
      return None
    if m < n - 1:
      return None
    if n * d < 2 * m:
      return None

    rng = random.Random(seed)

    pts = [(rng.random(), rng.random()) for _ in range(n)]
    points = [(x + EPS * rng.uniform(-1, 1), y + EPS * rng.uniform(-1, 1))
              for x, y in pts]
    tris = delaunay_bowyer_watson(points)
    all_edges = edges_from_tris(tris)

    rng.shuffle(all_edges)
    # edges_m = sorted(all_edges[:m])
    degrees = [0 for _ in range(n)]
    G = Graph(n)
    for u, v in all_edges:
      if len(G.edges) >= m:
        break
      if degrees[u] < d and degrees[v] < d:
        G.AddEdge(u, v)
        degrees[u] += 1
        degrees[v] += 1

    # for u, v in edges_m:
    #   G.AddEdge(u, v)
    return G

    # G = Graph(n)
    # C = int(math.ceil(math.sqrt(n)))
    # R = (n + C - 1) // C
    # lens = []
    # for r in range(R):
    #   left = n - r * C
    #   lens.append(math.min(C, math.max(0, left)))

    # def vid(r, c):
    #   if r < 0 or r >= R or c < 0 or c >= lens[r]:
    #     return -1
    #   idx = r * C + c
    #   return idx if idx < n else -1

    # cand = set()
    # def add(u, v):
    #   if u < 0 or v < 0 or u == v:
    #     return
    #   if u > v:
    #     u, v = v, u
    #   cand.add((u, v))

    # for r in range(R):
    #   for c in range(lens[r]):
    #     u = vid(r, c)
    #     if c + 1 < lens[r]:
    #       add(u, vid(r, c + 1))
    #     if r + 1 < R:
    #       if r & 1:
    #         if c < lens[r + 1]:
    #           add(u, vid(r + 1, c))
    #         if c + 1 < lens[r + 1]:
    #           add(u, vid(r + 1, c + 1))
    #       else:
    #         if c - 1 >= 0 and c - 1 < lens[r + 1]:
    #           add(u, vid(r + 1, c - 1))
    #         if c < lens[r + 1]:
    #           add(u, vid(r + 1, c))

    # cand = list(cand)
    # cnt = len(cand)
    # m_cap = math.min(d * n >> 1, cnt)
    # if m > m_cap:
    #   return None

    # deg = [0 for _ in range(n)]
    # def try_add(u, v, s):
    #   if (u, v) in s:
    #     return False
    #   if deg[u] >= d or deg[v] >= d:
    #     return False
    #   s.add((u, v))
    #   deg[u] += 1
    #   deg[v] += 1
    #   return True
