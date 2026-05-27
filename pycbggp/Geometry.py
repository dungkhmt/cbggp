import math

class Point:
  def __init__(self, x, y) -> None:
    self.x = x
    self.y = y

  def __lt__(self, p):
    return (self.x, self.y) < (p.x, p.y)
  def __eq__(self, p):
    return self.x == p.x and self.y == p.y
  
  def __add__(self, p):
    return Point(self.x + p.x, self.y + p.y)
  def __sub__(self, p):
    return Point(self.x - p.x, self.y - p.y)
  def __mul__(self, k):
    return Point(self.x * k, self.y * k)
  def __truediv__(self, k):
    return Point(self.x / k, self.y / k)
  
  def dot(self, p):
    return self.x * p.x + self.y * p.y
  def cross(self, p):
    return self.x * p.y - self.y * p.x
  def cross(self, p1, p2):
    return (p1 - self).cross(p2 - self)
  def dist2(self):
    return self.x * self.x + self.y * self.y
  def dist(self) -> float:
    return math.sqrt(self.dist2())
  
  def angle(self):
    return math.atan2(self.y, self.x)
  def unit(self):
    return self / self.dist()
  def perp(self):
    return Point(-self.y, self.x)
  def rotate(self, theta):
    c = math.cos(theta)
    s = math.sin(theta)
    return Point(self.x * c - self.y * s, self.x * s + self.y * c)

  def __repr__(self) -> str:
    return f"({self.x}, {self.y})"
