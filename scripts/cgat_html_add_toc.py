################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
'''
cgat_html_add_toc.py - insert table of contents in html document
================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Read an html document on stdin and add a table of contents
based on section headings in the document.

This document uses the ``<h1></h1>``, ``<h2></h2>``, ... html tags.

Usage
-----

Example::

   python cgat_html_add_toc.py --help

Type::

   python cgat_html_add_toc.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import re
import sys
import time
import getopt
import StringIO
from types import *

def usage():
    print """Usage: %s [--max-depth N] [--no-toc-links] file.html file-created.html
--max-depth N: The depth of the table of content should be N.
--no-toc-links: Don't insert links to the TOC after each heading.
--hide:
  Hide all headings and the content "below" which is marked with class="hide"
  Example: <h2 class="hide">Ignore Me</h2>
  
Number headings in a HTML file. The headings get numbered by
the <h?> tags in your file. You can insert a TOC (Table of Content)
with <!-- TOC -->.

  Example:
  <h1>First</h1>
  <h2>Sub-Item</h2>
  <h2>Sub-Item</h2>

  <h1>Second</h1>

  Result:

  <h1>1 First</h1>
  <h2>1.1 Sub-Item</h2>
  <h2>1.2 Sub-Item</h2>

  <h1>2 Second</h1>    """ % (
        os.path.basename(sys.argv[0]))

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["max-depth=",
                                                      "no-toc-links",
                                                      "hide",
                                                      "help",
                                                      "version"])
    except getopt.GetoptError, msg:
        usage()
        print msg
        sys.exit(1)

    max_depth=3
    no_toc_links=0
    do_hide=False
    for o, a in opts:
        if o == "--help":
            usage()
            sys.exit(0)
        elif o in ( "--version", ):
            print "version="
            sys.exit(0)
        elif o=="--no-toc-links":
            no_toc_links=1
        elif o=="--max-depth":
            try:
                max_depth=int(a)
            except ValueError:
                sys.stderr.write("%s must be a number\n" % a)
                sys.exit(1)
        elif o=="--hide":
            do_hide=True
        else:
            raise("Internal Error %s %s not processed" % (
                o, a))

    if len(args)!=2:
        usage()
        sys.exit(1)

    filename_in=args[0]
    filename_out=args[1]
    fd=open(filename_in)
    file=fd.readlines()
    fd.close()

    re_heading=re.compile(r'^(.*)<h(\d)([^>]*)>(.*?)<(/)?h(\d)>(.*)$')

    stack=[]
    last_level=0
    out=StringIO.StringIO()
    out_orig=out
    hide_out=StringIO.StringIO()
    hide=None
    headings=[]
    first = True
    for line in file:
        match=re_heading.match(line)
        if not match:
            out.write(line)
        else:
            heading=match.group(4)
            before_heading=match.group(1)
            attributes=match.group(3)
            
            after_heading=match.group(7)
            if match.group(2)!=match.group(6):
                raise "Parse Error: <h?> does not match </h?> %s %s: '%s' " % (
                    match.group(2), match.group(6), line)
            if match.group(5)!="/":
                raise "Parse Error: Missing slash in end-tag: '%s'" % line.strip()
            level=int(match.group(2))
            assert(level>0 and level<10)
            ## deal with TOCs, that do not start with <h1>
            if first:
                last_level = level - 1
                first = False
            if level==last_level+1:
                stack.append(1)
            elif level==last_level:
                stack[-1]+=1
            elif level<last_level:
                diff=last_level-level
                if diff>len(stack):
                    raise("Parse Error: diff: %s stack: %s last_level: %s"
                          " level: %s" % (
                        diff, stack, last_level, level))
                for i in range(diff):
                    stack.pop()
                stack[-1]+=1
            else:
                raise("Strange sequence of <h?> Stack: %s "
                      " level:%s last_level:%s line: '%s'" %
                      (stack, level, last_level, line))
            if hide and level<=hide:
                out=out_orig
                hide=None
            last_level=level
            number=[]
            space=[]
            number_td=[]
            for i in range(len(stack)):
                number.append(str(stack[i]))
                if i!=0:
                    space.append('&nbsp;&nbsp;')
            space=''.join(space)
            number='.'.join(number)

            match=re.match(r'^(.*?)\s*id="(.+?)"\s*(.*)$', attributes)
            if match:
                # If the heading has already an id, then take this name
                # and don't use link_N
                id=match.group(2)
                attributes="%s %s" % (match.group(1), match.group(3))
            else:
                id="link_%s" % number

            if do_hide:
                match=re.match(r'^(.*?)\s*class="hide"\s*(.*)$', attributes)
                if match:
                    # hide this heading and all stuff below
                    assert hide==None
                    hide=level
                    attributes="%s %s" % (match.group(1), match.group(2))
                    out=hide_out

            # Create column for each number in heading (for alignment)
            #print stack, heading, number
            if len(stack)==1:
                size="font-size: 100%"
            elif len(stack)==2:
                size="font-size: 80%"
            else:
                size="font-size: 60%"
            for i in range(max_depth):
                if i==len(stack)-1:
                    point="&nbsp;"
                    value=stack[i]
                elif i<len(stack):
                    point="."
                    value=stack[i]
                else:
                    point="&nbsp;"
                    value=""
                number_td.append(
                    '<td align="right" style="%s">%s%s</td>' % (
                    size, value, point))
            number_td=''.join(number_td)
            
            if no_toc_links:
                toc=""
            else:
                toc='<a style="font-size: 50%%" href="#toc_%s">[toc]</a>' % (
                    number)
                
            out.write(before_heading)
            out_orig.write('''
            <h%s %s id="%s">%s %s
             %s
            </h%s>''' % (level, attributes,
                   id, number, heading, toc,
                   level))
            if hide==level:
                out_orig.write("...")
            else:
                out.write(after_heading)
            assert(type(level)==IntType)

            #  create heading
            if level<=max_depth:
                headings.append('''
                 <tr>
                  %s
                  <td><a name="toc_%s" href="#%s">%s%s</a></td>
                 </tr>''' % (number_td, number,
                             id, space, heading))
                             
                
    headings='''
     <table>
      %s
     </table>''' % ''.join(headings)
    out=out.getvalue()
    out=re.sub(r'<!--\s*TOC\s*-->', ''.join(headings), out)
    out=re.sub(r'<html>',
               '''<html>
    <!--
     DO NOT EDIT!
    created by number-html-headings.py on %s from %s
     DO NOT EDIT!
     -->''' %
               (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), filename_in),
               out)
    out=re.sub(r'\s*</pre>', r'</pre>', out, re.DOTALL)
    fd=open(filename_out, "wt")
    fd.write(out)
    fd.close()

if __name__=="__main__":
    main()
