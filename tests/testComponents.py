import CGAT.Components as Components

c = Components.SComponents()

links = ( ( "1", "2"),
          ( "1", "3"),
          ( "2", "3"),
          ( "3", "3"),
          ( "4", "5"),
          ( "6", "6") )

for a,b in links:
    print a, b, c.add( a, b)

for x in "01234567":
    print x, c.get( x )

print c.getNumNodes()

print c.getComponents()

c = Components.IComponents()

for a,b in links:
    print a, b, c.add( int(a), int(b))

for x in range( 0, 8):
    print x, c.get( x )

print c.getNumNodes()

print c.getComponents()

c = Components.IComponents()
print c.getComponents()
print c.getNumNodes()

c = Components.IComponents()

for x in range( 3):
    c.add( x, x)
print c.getNumNodes()
print c.getComponents()



