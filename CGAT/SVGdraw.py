'''
SVGdraw.py - generate SVG drawings
======================================================

:Tags: Python

This module has been copied from 3rd party resources.

SVGdraw uses an object model drawing and a method toXML to create SVG graphics
by using easy to use classes and methods usualy you start by creating a drawing eg

    d=drawing()
    #then you create a SVG root element
    s=svg()
    #then you add some elements eg a circle and add it to the svg root element
    c=circle()
    #you can supply attributes by using named arguments.
    c=circle(fill='red',stroke='blue')
    #or by updating the attributes attribute:
    c.attributes['stroke-width']=1
    s.addElement(c)
    #then you add the svg root element to the drawing
    d.setSVG(s)
    #and finaly you xmlify the drawing
    d.toXml()
    

this results in the svg source of the drawing, which consists of a circle
on a white background. Its as easy as that;)
This module was created using the SVG specification of www.w3c.org and the
O'Reilly (www.oreilly.com) python books as information sources. A svg viewer
is available from www.adobe.com
'''

from six import StringIO
import sys

__version__ = "1.0"

use_dom_implementation = 0

if use_dom_implementation != 0:
    try:
        from xml.dom import implementation
        from xml.dom.ext import PrettyPrint
    except:
        raise ImportError("PyXML is required for using the dom implementation")


sys.setrecursionlimit = 50


def _escape(data, entities={}):
    """Escape &, <, and > in a string of data.

    You can escape other strings of data by passing a dictionary as
    the optional entities parameter.  The keys and values must all be
    strings; each key will be replaced with its corresponding value.
    """
    data = data.replace("&", "&amp;")
    data = data.replace("<", "&lt;")
    data = data.replace(">", "&gt;")
    for chars, entity in list(entities.items()):
        data = data.replace(chars, entity)
    return data


def _quoteattr(data, entities={}):
    """Escape and quote an attribute value.

    Escape &, <, and > in a string of data, then quote it for use as
    an attribute value.  The \" character will be escaped as well, if
    necessary.

    You can escape other strings of data by passing a dictionary as
    the optional entities parameter.  The keys and values must all be
    strings; each key will be replaced with its corresponding value.
    """
    data = _escape(data, entities)
    if '"' in data:
        if "'" in data:
            data = '"%s"' % data.replace('"', "&quot;")
        else:
            data = "'%s'" % data
    else:
        data = '"%s"' % data
    return data


def _xypointlist(a):
    """formats a list of xy pairs"""
    s = ''
    for e in a:  # this could be done more elegant
        s += str(e)[1:-1] + '  '
    return s


def _viewboxlist(a):
    """formats a tuple"""
    s = ''
    for e in a:
        s += str(e) + ' '
    return s


def _pointlist(a):
    """formats a list of numbers"""
    return str(a)[1:-1]


class pathdata:

    """class used to create a pathdata object which can be used for a path.
    although most methods are pretty straightforward it might be useful to look at the SVG specification."""
    # I didn't test the methods below.

    def __init__(self, x=None, y=None):
        self.path = []
        if x is not None and y is not None:
            self.path.append('M ' + str(x) + ' ' + str(y))

    def closepath(self):
        """ends the path"""
        self.path.append('z')

    def move(self, x, y):
        """move to absolute"""
        self.path.append('M ' + str(x) + ' ' + str(y))

    def relmove(self, x, y):
        """move to relative"""
        self.path.append('m ' + str(x) + ' ' + str(y))

    def line(self, x, y):
        """line to absolute"""
        self.path.append('L ' + str(x) + ' ' + str(y))

    def relline(self, x, y):
        """line to relative"""
        self.path.append('l ' + str(x) + ' ' + str(y))

    def hline(self, x):
        """horizontal line to absolute"""
        self.path.append('H' + str(x))

    def relhline(self, x):
        """horizontal line to relative"""
        self.path.append('h' + str(x))

    def vline(self, y):
        """verical line to absolute"""
        self.path.append('V' + str(y))

    def relvline(self, y):
        """vertical line to relative"""
        self.path.append('v' + str(y))

    def bezier(self, x1, y1, x2, y2, x, y):
        """bezier with xy1 and xy2 to xy absolut"""
        self.path.append('C' + str(x1) + ',' + str(y1) + ' ' +
                         str(x2) + ',' + str(y2) + ' ' + str(x) + ',' + str(y))

    def relbezier(self, x1, y1, x2, y2, x, y):
        """bezier with xy1 and xy2 to xy relative"""
        self.path.append('c' + str(x1) + ',' + str(y1) + ' ' +
                         str(x2) + ',' + str(y2) + ' ' + str(x) + ',' + str(y))

    def smbezier(self, x2, y2, x, y):
        """smooth bezier with xy2 to xy absolut"""
        self.path.append(
            'S' + str(x2) + ',' + str(y2) + ' ' + str(x) + ',' + str(y))

    def relsmbezier(self, x2, y2, x, y):
        """smooth bezier with xy2 to xy relative"""
        self.path.append(
            's' + str(x2) + ',' + str(y2) + ' ' + str(x) + ',' + str(y))

    def qbezier(self, x1, y1, x, y):
        """quadratic bezier with xy1 to xy absolut"""
        self.path.append(
            'Q' + str(x1) + ',' + str(y1) + ' ' + str(x) + ',' + str(y))

    def relqbezier(self, x1, y1, x, y):
        """quadratic bezier with xy1 to xy relative"""
        self.path.append(
            'q' + str(x1) + ',' + str(y1) + ' ' + str(x) + ',' + str(y))

    def smqbezier(self, x, y):
        """smooth quadratic bezier to xy absolut"""
        self.path.append('T' + str(x) + ',' + str(y))

    def relsmqbezier(self, x, y):
        """smooth quadratic bezier to xy relative"""
        self.path.append('t' + str(x) + ',' + str(y))

    def ellarc(self, rx, ry, xrot, laf, sf, x, y):
        """elliptival arc with rx and ry rotating with xrot using large-arc-flag and sweep-flag  to xy absolut"""
        self.path.append('A' + str(rx) + ',' + str(ry) + ' ' + str(xrot) +
                         ' ' + str(laf) + ' ' + str(sf) + ' ' + str(x) + ' ' + str(y))

    def relellarc(self, rx, ry, xrot, laf, sf, x, y):
        """elliptival arc with rx and ry rotating with xrot using large-arc-flag and sweep-flag  to xy relative"""
        self.path.append('a' + str(rx) + ',' + str(ry) + ' ' + str(xrot) +
                         ' ' + str(laf) + ' ' + str(sf) + ' ' + str(x) + ' ' + str(y))

    def __repr__(self):
        return ' '.join(self.path)


class SVGelement:

    """SVGelement(type,attributes,elements,text,namespace,**args)
    Creates a arbitrary svg element and is intended to be subclassed not used on its own.
    This element is the base of every svg element it defines a class which resembles
    a xml-element. The main advantage of this kind of implementation is that you don't
    have to create a toXML method for every different graph object. Every element
    consists of a type, attribute, optional subelements, optional text and an optional
    namespace. Note the elements==None, if elements = None:self.elements=[] construction.
    This is done because if you default to elements=[] every object has a reference
    to the same empty list."""

    def __init__(self, type='', attributes=None, elements=None, text='', namespace='', cdata=None, **args):
        self.type = type
        if attributes is None:
            self.attributes = {}
        else:
            self.attributes = attributes
        if elements is None:
            self.elements = []
        else:
            self.elements = elements
        self.text = text
        self.namespace = namespace
        self.cdata = cdata
        for arg in list(args.keys()):
            self.attributes[arg] = args[arg]

    def addElement(self, SVGelement):
        """adds an element to a SVGelement

        SVGelement.addElement(SVGelement)
        """
        self.elements.append(SVGelement)

    def toXml(self, level, f):
        f.write('\t' * level)
        f.write('<' + self.type)
        for attkey in list(self.attributes.keys()):
            f.write(' ' + _escape(str(attkey)) + '=' +
                    _quoteattr(str(self.attributes[attkey])))
        if self.namespace:
            f.write(' xmlns="' + _escape(str(self.namespace)) + '" ')
            # added by Andreas Heger
            f.write(
                ' xmlns:xlink="' + _escape(str("http://www.w3.org/1999/xlink")) + '" ')
        if self.elements or self.text or self.cdata:
            f.write('>')
        if self.elements:
            f.write('\n')
        for element in self.elements:
            element.toXml(level + 1, f)
        if self.cdata:
            f.write('\n' + '\t' * (level + 1) + '<![CDATA[')
            for line in self.cdata.splitlines():
                f.write('\n' + '\t' * (level + 2) + line)
            f.write('\n' + '\t' * (level + 1) + ']]>\n')
        if self.text:
            if isinstance(self.text, str):  # If the text is only text
                f.write(_escape(str(self.text)))
            else:  # If the text is a spannedtext class
                f.write(str(self.text))
        if self.elements:
            f.write('\t' * level + '</' + self.type + '>\n')
        elif self.text:
            f.write('</' + self.type + '>\n')
        elif self.cdata:
            f.write('\t' * level + '</' + self.type + '>\n')
        else:
            f.write('/>\n')


class tspan(SVGelement):

    """ts=tspan(text='',**args)

    a tspan element can be used for applying formatting to a textsection
    usage:
    ts=tspan('this text is bold')
    ts.attributes['font-weight']='bold'
    st=spannedtext()
    st.addtspan(ts)
    t=text(3,5,st)
    """

    def __init__(self, text=None, **args):
        SVGelement.__init__(self, 'tspan', **args)
        if self.text is not None:
            self.text = text

    def __repr__(self):
        s = "<tspan"
        for key, value in list(self.attributes.items()):
            s += ' %s="%s"' % (key, value)
        s += '>'
        s += self.text
        s += '</tspan>'
        return s


class tref(SVGelement):

    """tr=tref(link='',**args)

    a tref element can be used for referencing text by a link to its id.
    usage:
    tr=tref('#linktotext')
    st=spannedtext()
    st.addtref(tr)
    t=text(3,5,st)
    """

    def __init__(self, link, **args):
        SVGelement.__init__(self, 'tref', {'xlink:href': link}, **args)

    def __repr__(self):
        s = "<tref"

        for key, value in list(self.attributes.items()):
            s += ' %s="%s"' % (key, value)
        s += '/>'
        return s


class spannedtext:

    """st=spannedtext(textlist=[])

    a spannedtext can be used for text which consists of text, tspan's and tref's
    You can use it to add to a text element or path element. Don't add it directly
    to a svg or a group element.
    usage:

    ts=tspan('this text is bold')
    ts.attributes['font-weight']='bold'
    tr=tref('#linktotext')
    tr.attributes['fill']='red'
    st=spannedtext()
    st.addtspan(ts)
    st.addtref(tr)
    st.addtext('This text is not bold')
    t=text(3,5,st)
    """

    def __init__(self, textlist=None):
        if textlist is None:
            self.textlist = []
        else:
            self.textlist = textlist

    def addtext(self, text=''):
        self.textlist.append(text)

    def addtspan(self, tspan):
        self.textlist.append(tspan)

    def addtref(self, tref):
        self.textlist.append(tref)

    def __repr__(self):
        s = ""
        for element in self.textlist:
            s += str(element)
        return s


class rect(SVGelement):

    """r=rect(width,height,x,y,fill,stroke,stroke_width,**args)

    a rectangle is defined by a width and height and a xy pair 
    """

    def __init__(self, x=None, y=None, width=None, height=None, fill=None, stroke=None, stroke_width=None, **args):
        if width is None or height is None:
            if width is not None:
                raise ValueError('height is required')
            if height is not None:
                raise ValueError('width is required')
            else:
                raise ValueError('both height and width are required')
        SVGelement.__init__(
            self, 'rect', {'width': width, 'height': height}, **args)
        if x is not None:
            self.attributes['x'] = x
        if y is not None:
            self.attributes['y'] = y
        if fill is not None:
            self.attributes['fill'] = fill
        if stroke is not None:
            self.attributes['stroke'] = stroke
        if stroke_width is not None:
            self.attributes['stroke-width'] = stroke_width


class ellipse(SVGelement):

    """e=ellipse(rx,ry,x,y,fill,stroke,stroke_width,**args)

    an ellipse is defined as a center and a x and y radius.
    """

    def __init__(self, cx=None, cy=None, rx=None, ry=None, fill=None, stroke=None, stroke_width=None, **args):
        if rx is None or ry is None:
            if rx is not None:
                raise ValueError('rx is required')
            if ry is not None:
                raise ValueError('ry is required')
            else:
                raise ValueError('both rx and ry are required')
        SVGelement.__init__(self, 'ellipse', {'rx': rx, 'ry': ry}, **args)
        if cx is not None:
            self.attributes['cx'] = cx
        if cy is not None:
            self.attributes['cy'] = cy
        if fill is not None:
            self.attributes['fill'] = fill
        if stroke is not None:
            self.attributes['stroke'] = stroke
        if stroke_width is not None:
            self.attributes['stroke-width'] = stroke_width


class circle(SVGelement):

    """c=circle(x,y,radius,fill,stroke,stroke_width,**args)

    The circle creates an element using a x, y and radius values eg
    """

    def __init__(self, cx=None, cy=None, r=None, fill=None, stroke=None, stroke_width=None, **args):
        if r is None:
            raise ValueError('r is required')
        SVGelement.__init__(self, 'circle', {'r': r}, **args)
        if cx is not None:
            self.attributes['cx'] = cx
        if cy is not None:
            self.attributes['cy'] = cy
        if fill is not None:
            self.attributes['fill'] = fill
        if stroke is not None:
            self.attributes['stroke'] = stroke
        if stroke_width is not None:
            self.attributes['stroke-width'] = stroke_width


class point(circle):

    """p=point(x,y,color)

    A point is defined as a circle with a size 1 radius. It may be more efficient to use a
    very small rectangle if you use many points because a circle is difficult to render.
    """

    def __init__(self, x, y, fill='black', **args):
        circle.__init__(self, x, y, 1, fill, **args)


class line(SVGelement):

    """l=line(x1,y1,x2,y2,stroke,stroke_width,**args)

    A line is defined by a begin x,y pair and an end x,y pair
    """

    def __init__(self, x1=None, y1=None, x2=None, y2=None, stroke=None, stroke_width=None, **args):
        SVGelement.__init__(self, 'line', **args)
        if x1 is not None:
            self.attributes['x1'] = x1
        if y1 is not None:
            self.attributes['y1'] = y1
        if x2 is not None:
            self.attributes['x2'] = x2
        if y2 is not None:
            self.attributes['y2'] = y2
        if stroke_width is not None:
            self.attributes['stroke-width'] = stroke_width
        if stroke is not None:
            self.attributes['stroke'] = stroke


class polyline(SVGelement):

    """pl=polyline([[x1,y1],[x2,y2],...],fill,stroke,stroke_width,**args)

    a polyline is defined by a list of xy pairs
    """

    def __init__(self, points, fill=None, stroke=None, stroke_width=None, **args):
        SVGelement.__init__(
            self, 'polyline', {'points': _xypointlist(points)}, **args)
        if fill is not None:
            self.attributes['fill'] = fill
        if stroke_width is not None:
            self.attributes['stroke-width'] = stroke_width
        if stroke is not None:
            self.attributes['stroke'] = stroke


class polygon(SVGelement):

    """pl=polyline([[x1,y1],[x2,y2],...],fill,stroke,stroke_width,**args)

    a polygon is defined by a list of xy pairs
    """

    def __init__(self, points, fill=None, stroke=None, stroke_width=None, **args):
        SVGelement.__init__(
            self, 'polygon', {'points': _xypointlist(points)}, **args)
        if fill is not None:
            self.attributes['fill'] = fill
        if stroke_width is not None:
            self.attributes['stroke-width'] = stroke_width
        if stroke is not None:
            self.attributes['stroke'] = stroke


class path(SVGelement):

    """p=path(path,fill,stroke,stroke_width,**args)

    a path is defined by a path object and optional width, stroke and fillcolor
    """

    def __init__(self, pathdata, fill=None, stroke=None, stroke_width=None, id=None, **args):
        SVGelement.__init__(self, 'path', {'d': str(pathdata)}, **args)
        if stroke is not None:
            self.attributes['stroke'] = stroke
        if fill is not None:
            self.attributes['fill'] = fill
        if stroke_width is not None:
            self.attributes['stroke-width'] = stroke_width
        if id is not None:
            self.attributes['id'] = id


class text(SVGelement):

    """t=text(x,y,text,font_size,font_family,**args)

    a text element can bge used for displaying text on the screen
    """

    def __init__(self, x=None, y=None, text=None, font_size=None, font_family=None, text_anchor=None, font_style=None, **args):
        SVGelement.__init__(self, 'text', **args)
        if x is not None:
            self.attributes['x'] = x
        if y is not None:
            self.attributes['y'] = y
        if font_size is not None:
            self.attributes['font-size'] = font_size
        if font_family is not None:
            self.attributes['font-family'] = font_family
        if font_style is not None:
            self.attributes['font-style'] = font_style
        if text is not None:
            self.text = text
        if text_anchor is not None:
            self.attributes['text-anchor'] = text_anchor


class textpath(SVGelement):

    """tp=textpath(text,link,**args)

    a textpath places a text on a path which is referenced by a link.   
    """

    def __init__(self, link, text=None, **args):
        SVGelement.__init__(self, 'textPath', {'xlink:href': link}, **args)
        if text is not None:
            self.text = text


class pattern(SVGelement):

    """p=pattern(x,y,width,height,patternUnits,**args)

    A pattern is used to fill or stroke an object using a pre-defined
    graphic object which can be replicated ("tiled") at fixed intervals
    in x and y to cover the areas to be painted.
    """

    def __init__(self, x=None, y=None, width=None, height=None, patternUnits=None, **args):
        SVGelement.__init__(self, 'pattern', **args)
        if x is not None:
            self.attributes['x'] = x
        if y is not None:
            self.attributes['y'] = y
        if width is not None:
            self.attributes['width'] = width
        if height is not None:
            self.attributes['height'] = height
        if patternUnits is not None:
            self.attributes['patternUnits'] = patternUnits


class title(SVGelement):

    """t=title(text,**args)

    a title is a text element. The text is displayed in the title bar
    add at least one to the root svg element
    """

    def __init__(self, text=None, **args):
        SVGelement.__init__(self, 'title', **args)
        if text is not None:
            self.text = text


class description(SVGelement):

    """d=description(text,**args)

    a description can be added to any element and is used for a tooltip
    Add this element before adding other elements.
    """

    def __init__(self, text=None, **args):
        SVGelement.__init__(self, 'desc', **args)
        if text is not None:
            self.text = text


class lineargradient(SVGelement):

    """lg=lineargradient(x1,y1,x2,y2,id,**args)

    defines a lineargradient using two xy pairs.
    stop elements van be added to define the gradient colors.
    """

    def __init__(self, x1=None, y1=None, x2=None, y2=None, id=None, **args):
        SVGelement.__init__(self, 'linearGradient', **args)
        if x1 is not None:
            self.attributes['x1'] = x1
        if y1 is not None:
            self.attributes['y1'] = y1
        if x2 is not None:
            self.attributes['x2'] = x2
        if y2 is not None:
            self.attributes['y2'] = y2
        if id is not None:
            self.attributes['id'] = id


class radialgradient(SVGelement):

    """rg=radialgradient(cx,cy,r,fx,fy,id,**args)

    defines a radial gradient using a outer circle which are defined by a cx,cy and r and by using a focalpoint.
    stop elements van be added to define the gradient colors.
    """

    def __init__(self, cx=None, cy=None, r=None, fx=None, fy=None, id=None, **args):
        SVGelement.__init__(self, 'radialGradient', **args)
        if cx is not None:
            self.attributes['cx'] = cx
        if cy is not None:
            self.attributes['cy'] = cy
        if r is not None:
            self.attributes['r'] = r
        if fx is not None:
            self.attributes['fx'] = fx
        if fy is not None:
            self.attributes['fy'] = fy
        if id is not None:
            self.attributes['id'] = id


class stop(SVGelement):

    """st=stop(offset,stop_color,**args)

    Puts a stop color at the specified radius
    """

    def __init__(self, offset, stop_color=None, **args):
        SVGelement.__init__(self, 'stop', {'offset': offset}, **args)
        if stop_color is not None:
            self.attributes['stop-color'] = stop_color


class style(SVGelement):

    """st=style(type,cdata=None,**args)

    Add a CDATA element to this element for defing in line stylesheets etc..
    """

    def __init__(self, type, cdata=None, **args):
        SVGelement.__init__(self, 'style', {'type': type}, cdata=cdata, **args)


class image(SVGelement):

    """im=image(url,width,height,x,y,**args)

    adds an image to the drawing. Supported formats are .png, .jpg and .svg.
    """

    def __init__(self, url, x=None, y=None, width=None, height=None, **args):
        if width is None or height is None:
            if width is not None:
                raise ValueError('height is required')
            if height is not None:
                raise ValueError('width is required')
            else:
                raise ValueError('both height and width are required')
        SVGelement.__init__(
            self, 'image', {'xlink:href': url, 'width': width, 'height': height}, **args)
        if x is not None:
            self.attributes['x'] = x
        if y is not None:
            self.attributes['y'] = y


class cursor(SVGelement):

    """c=cursor(url,**args)

    defines a custom cursor for a element or a drawing
    """

    def __init__(self, url, **args):
        SVGelement.__init__(self, 'cursor', {'xlink:href': url}, **args)


class marker(SVGelement):

    """m=marker(id,viewbox,refX,refY,markerWidth,markerHeight,**args)

    defines a marker which can be used as an endpoint for a line or other pathtypes
    add an element to it which should be used as a marker.
    """

    def __init__(self, id=None, viewBox=None, refx=None, refy=None, markerWidth=None, markerHeight=None, **args):
        SVGelement.__init__(self, 'marker', **args)
        if id is not None:
            self.attributes['id'] = id
        if viewBox is not None:
            self.attributes['viewBox'] = _viewboxlist(viewBox)
        if refx is not None:
            self.attributes['refX'] = refx
        if refy is not None:
            self.attributes['refY'] = refy
        if markerWidth is not None:
            self.attributes['markerWidth'] = markerWidth
        if markerHeight is not None:
            self.attributes['markerHeight'] = markerHeight


class group(SVGelement):

    """g=group(id,**args)

    a group is defined by an id and is used to contain elements
    g.addElement(SVGelement)
    """

    def __init__(self, id=None, **args):
        SVGelement.__init__(self, 'g', **args)
        if id is not None:
            self.attributes['id'] = id


class symbol(SVGelement):

    """sy=symbol(id,viewbox,**args)

    defines a symbol which can be used on different places in your graph using
    the use element. A symbol is not rendered but you can use 'use' elements to
    display it by referencing its id.
    sy.addElement(SVGelement)
    """

    def __init__(self, id=None, viewBox=None, **args):
        SVGelement.__init__(self, 'symbol', **args)
        if id is not None:
            self.attributes['id'] = id
        if viewBox is not None:
            self.attributes['viewBox'] = _viewboxlist(viewBox)


class defs(SVGelement):

    """d=defs(``**args``)

    container for defining elements
    """

    def __init__(self, **args):
        SVGelement.__init__(self, 'defs', **args)


class switch(SVGelement):

    """sw=switch(``**args``)

    Elements added to a switch element which are "switched" by the attributes
    requiredFeatures, requiredExtensions and systemLanguage.
    Refer to the SVG specification for details.
    """

    def __init__(self, **args):
        SVGelement.__init__(self, 'switch', **args)


class use(SVGelement):

    """u=use(link,x,y,width,height,``**args``)

    references a symbol by linking to its id and its position, height and width
    """

    def __init__(self, link, x=None, y=None, width=None, height=None, **args):
        SVGelement.__init__(self, 'use', {'xlink:href': link}, **args)
        if x is not None:
            self.attributes['x'] = x
        if y is not None:
            self.attributes['y'] = y

        if width is not None:
            self.attributes['width'] = width
        if height is not None:
            self.attributes['height'] = height


class link(SVGelement):

    """a=link(url,``**args``)

    a link  is defined by a hyperlink. add elements which have to be linked
    a.addElement(SVGelement)
    """

    def __init__(self, link='', **args):
        SVGelement.__init__(self, 'a', {'xlink:href': link}, **args)


class view(SVGelement):

    """v=view(id,``**args``)

    a view can be used to create a view with different attributes"""

    def __init__(self, id=None, **args):
        SVGelement.__init__(self, 'view', **args)
        if id is not None:
            self.attributes['id'] = id


class script(SVGelement):

    """sc=script(type,type,cdata,``**args``)

    adds a script element which contains CDATA to the SVG drawing

    """

    def __init__(self, type, cdata=None, **args):
        SVGelement.__init__(
            self, 'script', {'type': type}, cdata=cdata, **args)


class animate(SVGelement):

    """an=animate(attribute,from,to,during,``**args``)

    animates an attribute.    
    """

    def __init__(self, attribute, fr=None, to=None, dur=None, **args):
        SVGelement.__init__(
            self, 'animate', {'attributeName': attribute}, **args)
        if fr is not None:
            self.attributes['from'] = fr
        if to is not None:
            self.attributes['to'] = to
        if dur is not None:
            self.attributes['dur'] = dur


class animateMotion(SVGelement):

    """an=animateMotion(pathdata,dur,``**args``)

    animates a SVGelement over the given path in dur seconds
    """

    def __init__(self, pathdata, dur, **args):
        SVGelement.__init__(self, 'animateMotion', **args)
        if pathdata is not None:
            self.attributes['path'] = str(pathdata)
        if dur is not None:
            self.attributes['dur'] = dur


class animateTransform(SVGelement):

    """antr=animateTransform(type,from,to,dur,``**args``)

    transform an element from and to a value.
    """

    def __init__(self, type=None, fr=None, to=None, dur=None, **args):
        SVGelement.__init__(
            self, 'animateTransform', {'attributeName': 'transform'}, **args)
        # As far as I know the attributeName is always transform
        if type is not None:
            self.attributes['type'] = type
        if fr is not None:
            self.attributes['from'] = fr
        if to is not None:
            self.attributes['to'] = to
        if dur is not None:
            self.attributes['dur'] = dur


class animateColor(SVGelement):

    """ac=animateColor(attribute,type,from,to,dur,``**args``)

    Animates the color of a element
    """

    def __init__(self, attribute, type=None, fr=None, to=None, dur=None, **args):
        SVGelement.__init__(
            self, 'animateColor', {'attributeName': attribute}, **args)
        if type is not None:
            self.attributes['type'] = type
        if fr is not None:
            self.attributes['from'] = fr
        if to is not None:
            self.attributes['to'] = to
        if dur is not None:
            self.attributes['dur'] = dur


class set(SVGelement):

    """st=set(attribute,to,during,``**args``)

    sets an attribute to a value for a
    """

    def __init__(self, attribute, to=None, dur=None, **args):
        SVGelement.__init__(self, 'set', {'attributeName': attribute}, **args)
        if to is not None:
            self.attributes['to'] = to
        if dur is not None:
            self.attributes['dur'] = dur


class svg(SVGelement):

    """s=svg(viewbox,width,height,``**args``)

    a svg or element is the root of a drawing add all elements to a svg element.
    You can have different svg elements in one svg file
    s.addElement(SVGelement)

    eg
    d=drawing()
    s=svg((0,0,100,100),'100%','100%')
    c=circle(50,50,20)
    s.addElement(c)
    d.setSVG(s)
    d.toXml()
    """

    def __init__(self, viewBox=None, width=None, height=None, **args):
        SVGelement.__init__(self, 'svg', **args)
        if viewBox is not None:
            self.attributes['viewBox'] = _viewboxlist(viewBox)
        if width is not None:
            self.attributes['width'] = width
        if height is not None:
            self.attributes['height'] = height
        self.namespace = "http://www.w3.org/2000/svg"


class drawing:

    """d=drawing()

    this is the actual SVG document. It needs a svg element as a root.
    Use the addSVG method to set the svg to the root. Use the toXml method to write the SVG
    source to the screen or to a file
    d=drawing()
    d.addSVG(svg)
    d.toXml(optionalfilename)
    """

    def __init__(self):
        self.svg = None

    def setSVG(self, svg):
        self.svg = svg
        # Voeg een element toe aan de grafiek toe.
    if use_dom_implementation == 0:
        def toXml(self, filename='', compress=False):
            xml = StringIO()
            xml.write("<?xml version='1.0' encoding='UTF-8'?>\n")
            xml.write(
                "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\" \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd \">\n")
            self.svg.toXml(0, xml)
            if not filename:
                if compress:
                    import gzip
                    f = StringIO()
                    zf = gzip.GzipFile(fileobj=f, mode='wb')
                    zf.write(xml.getvalue())
                    zf.close()
                    f.seek(0)
                    return f.read()
                else:
                    return xml.getvalue()
            else:
                if filename[-4:] == 'svgz':
                    import gzip
                    f = gzip.GzipFile(
                        filename=filename, mode="wb", compresslevel=9)
                    f.write(xml.getvalue())
                    f.close()
                else:
                    f = file(filename, 'w')
                    f.write(xml.getvalue())
                    f.close()

    else:
        def toXml(self, filename='', compress=False):
            """drawing.toXml()        ---->to the screen
            drawing.toXml(filename)---->to the file
            writes a svg drawing to the screen or to a file
            compresses if filename ends with svgz or if compress is true
            """
            doctype = implementation.createDocumentType(
                'svg', "-//W3C//DTD SVG 1.0//EN""", 'http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd ')

            global root
            # root is defined global so it can be used by the appender. Its also possible to use it as an arugument but
            # that is a bit messy.
            root = implementation.createDocument(None, None, doctype)
            # Create the xml document.
            global appender

            def appender(element, elementroot):
                """This recursive function appends elements to an element and sets the attributes
                and type. It stops when alle elements have been appended"""
                if element.namespace:
                    e = root.createElementNS(element.namespace, element.type)
                else:
                    e = root.createElement(element.type)
                if element.text:
                    textnode = root.createTextNode(element.text)
                    e.appendChild(textnode)
                # in element.attributes is supported from python 2.2
                for attribute in list(element.attributes.keys()):
                    e.setAttribute(
                        attribute, str(element.attributes[attribute]))
                if element.elements:
                    for el in element.elements:
                        e = appender(el, e)
                elementroot.appendChild(e)
                return elementroot
            root = appender(self.svg, root)
            if not filename:
                xml = StringIO()
                PrettyPrint(root, xml)
                if compress:
                    import gzip
                    f = StringIO()
                    zf = gzip.GzipFile(fileobj=f, mode='wb')
                    zf.write(xml.getvalue())
                    zf.close()
                    f.seek(0)
                    return f.read()
                else:
                    return xml.getvalue()
            else:
                try:
                    if filename[-4:] == 'svgz':
                        import gzip
                        xml = StringIO()
                        PrettyPrint(root, xml)
                        f = gzip.GzipFile(
                            filename=filename, mode='wb', compresslevel=9)
                        f.write(xml.getvalue())
                        f.close()
                    else:
                        f = open(filename, 'w')
                        PrettyPrint(root, f)
                        f.close()
                except:
                    print("Cannot write SVG file: " + filename)

    def validate(self):
        try:
            import xml.parsers.xmlproc.xmlval
        except:
            raise ImportError('PyXml is required for validating SVG')
        svg = self.toXml()
        xv = xml.parsers.xmlproc.xmlval.XMLValidator()
        try:
            xv.feed(svg)
        except:
            raise ValueError("SVG is not well formed, see messages above")
        else:
            print("SVG well formed")


if __name__ == '__main__':

    d = drawing()
    s = svg((0, 0, 100, 100))
    r = rect(-100, -100, 300, 300, 'cyan')
    s.addElement(r)

    t = title('SVGdraw Demo')
    s.addElement(t)
    g = group('animations')
    e = ellipse(0, 0, 5, 2)
    g.addElement(e)
    c = circle(0, 0, 1, 'red')
    g.addElement(c)
    pd = pathdata(0, -10)
    for i in range(6):
        pd.relsmbezier(10, 5, 0, 10)
        pd.relsmbezier(-10, 5, 0, 10)
    an = animateMotion(pd, 10)
    an.attributes['rotate'] = 'auto-reverse'
    an.attributes['repeatCount'] = "indefinite"
    g.addElement(an)
    s.addElement(g)
    for i in range(20, 120, 20):
        u = use('#animations', i, 0)
        s.addElement(u)
    for i in range(0, 120, 20):
        for j in range(5, 105, 10):
            c = circle(i, j, 1, 'red', 'black', .5)
            s.addElement(c)
    d.setSVG(s)

    print(d.toXml())
