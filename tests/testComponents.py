import CGAT.Components as Components

c = Components.SComponents()

# Todo: change Components library to enable unicode strings
links = ((b"1", b"2"),
         (b"1", b"3"),
         (b"2", b"3"),
         (b"3", b"3"),
         (b"4", b"5"),
         (b"6", b"6"))

for a, b in links:
    print(a, b, c.add(a, b))

for x in b"01234567":
    print(x, c.get(bytes(x)))

print(c.getNumNodes())

print(c.getComponents())

c = Components.IComponents()

for a, b in links:
    print(a, b, c.add(int(a), int(b)))

for x in range(0, 8):
    print(x, c.get(x))

print(c.getNumNodes())

print(c.getComponents())

c = Components.IComponents()
print(c.getComponents())
print(c.getNumNodes())

c = Components.IComponents()

for x in range(3):
    c.add(x, x)
print(c.getNumNodes())
print(c.getComponents())
