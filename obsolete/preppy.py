#copyright ReportLab Inc. 2000-2006
#see license.txt for license details

"""preppy - a Python preprocessor.

This is the Python equivalent of ASP or JSP - a preprocessor which lets you
embed python expressions, loops and conditionals, and 'scriptlets' in any
kind of text file.  It provides a very natural solution for generating
dynamic HTML pages, which is not connected to any particular web server
architecture.

You create a template file (conventionally ending in .prep) containing
python expressions, loops and conditionals, and scripts.  These occur
between double curly braces:
   Dear {{surname}},
   You owe us {{amount}} {{if amount>1000}}which is pretty serious{{endif}}

On first use or any any change in the template, this is normally converted to a
python source module 'in memory', then to a compiled pyc file which is saved to
disk alongside the original.  Options control this; you can operate entirely
in memory, or look at the generated python code if you wish.

On subsequent use, the generated module is imported and loaded directly.
The module contains a run(...) function; you can pass in a dictionary of
parameters (such as the surname and amount parameters above), and optionally
an output stream or output-collection function if you don't want it to go to
standard output.

The command line options let you run modules with hand-input parameters -
useful for basic testing - and also to batch-compile or clean directories.
As with python scripts, it is a good idea to compile prep files on installation,
since unix applications may run as a different user and not have the needed
permission to store compiled modules.

"""
VERSION = 0.8
__version__=''' $Id: preppy.py 2782 2009-09-10 11:40:29Z andreas $ '''

USAGE = """
The command line interface lets you test, compile and clean up:

    preppy modulename [arg1=value1, arg2=value2.....]
       - shorthand for 'preppy run ...', see below.

    preppy run modulename [arg1=value1, arg2=value2.....]
       - runs the module, optionally with arguments.  e.g.
         preppy.py flintstone.prep name=fred sex=m

    preppy.py compile [-f] [-v] [-p] module1[.prep] module2[.prep] module3 ...
       - compiles explicit modules

    preppy.py compile [-f] [-v] [-p] dirname1  dirname2 ...
       - compiles all prep files in directory recursively

    preppy.py clean dirname1 dirname2 ...19
       - removes any py or pyc files created from past compilations


"""

STARTDELIMITER = "{{"
ENDDELIMITER = "}}"
QSTARTDELIMITER = "{${"
QENDDELIMITER = "}$}"
QUOTE = "$"
QUOTEQUOTE = "$$"
# SEQUENCE OF REPLACEMENTS FOR UNESCAPING A STRING.
UNESCAPES = ((QSTARTDELIMITER, STARTDELIMITER), (QENDDELIMITER, ENDDELIMITER), (QUOTEQUOTE, QUOTE))

import re
import sys
import os
import imp
import struct
try:
    from hashlib import md5
except ImportError:
    from md5 import md5
from compiler import pycodegen, pyassem, future, consts

class _preppy_FunctionCodeGenerator(pycodegen.FunctionCodeGenerator):
    def __init__(self, func, scopes, isLambda, class_name, mod):
        self.scopes = scopes
        self.scope = scopes[func]
        self.class_name = class_name
        self.module = mod
        if isLambda:
            klass = FunctionCodeGenerator
            name = "<lambda.%d>" % klass.lambdaCount
            klass.lambdaCount = klass.lambdaCount + 1
        else:
            name = func.name
        args, hasTupleArg = pycodegen.generateArgList(func.argnames)
        self.optimized = int(name!='__code__')
        self._special_locals = name=='__code__'
        self.graph = pyassem.PyFlowGraph(name, func.filename, args,
                                         optimized=self.optimized)
        if name=='__code__': self.graph.setFlag(consts.CO_NEWLOCALS)
        self.isLambda = isLambda
        self.super_init()

        if not isLambda and func.doc:
            self.setDocstring(func.doc)

        lnf = pycodegen.walk(func.code, self.NameFinder(args), verbose=0)
        self.locals.push(lnf.getLocals())
        if func.varargs:
            self.graph.setFlag(consts.CO_VARARGS)
        if func.kwargs:
            self.graph.setFlag(consts.CO_VARKEYWORDS)
        self.set_lineno(func)
        if hasTupleArg:
            self.generateArgUnpack(func.argnames)
        self.graph.setFreeVars(self.scope.get_free_vars())
        self.graph.setCellVars(self.scope.get_cell_vars())
        if self.scope.generator is not None:
            self.graph.setFlag(pycodegen.CO_GENERATOR)

    def _nameOp(self, prefix, name):
        name = self.mangle(name)
        scope = self.scope.check_name(name)
        if self._special_locals and name=='locals': name = 'globals'
        if scope == consts.SC_LOCAL:
            if not (self.optimized or name in ['dictionary','__write__','__swrite__','outputfile', '__save_sys_stdout__', '__d__']):
                self.emit(prefix + '_NAME', name)
            else:
                self.emit(prefix + '_FAST', name)
        elif scope == consts.SC_GLOBAL:
            if not self.optimized:
                self.emit(prefix + '_NAME', name)
            else:
                self.emit(prefix + '_GLOBAL', name)
        elif scope == consts.SC_FREE or scope == consts.SC_CELL:
            self.emit(prefix + '_DEREF', name)
        else:
            raise RuntimeError, "unsupported scope for var %s: %d" % \
                  (name, scope)

class _preppy_ModuleCodeGenerator(pycodegen.ModuleCodeGenerator):
    def __init__(self, tree):
        self.graph = pyassem.PyFlowGraph("<module>", tree.filename)
        self.futures = future.find_futures(tree)
        pycodegen.ModuleCodeGenerator._ModuleCodeGenerator__super_init(self)
        self.__class__.FunctionGen = _preppy_FunctionCodeGenerator
        pycodegen.walk(tree, self)

def unescape(s, unescapes=UNESCAPES):
    for (old, new) in unescapes:
        s = s.replace(old, new)
    return s


"""
{{if x:}} this text $(endif}} that

becomes

if x:
    print (" this text"),
print (that\n")

note that extra spaces may be introduced compliments of print.
could change this to a file.write for more exact control...
"""

teststring = """
this test script should produce a runnable program
{{script}}
  class X:
      pass
  x = X()
  x.a = "THE A VALUE OF X"
  yislonger = "y is longer!"
  import math
  a = dictionary = {"key": "value", "key2": "value2", "10%": "TEN PERCENT"}
  loop = "LOOP"
{{endscript}}
this line has a percent in it 10%
here is the a value in x: {{x.a}}
just a norml value here: {{yislonger}} string {{a["10%"]}}
 the sine of 12.3 is {{math.sin(12.3)}}
 {{script}} a=0 {{endscript}}
 these parens should be empty
 ({{if a:}}
conditional text{{endif}})
 {{script}} a=1
 {{endscript}}
 these parens should be full
 ({{if a:}}
conditional text{{endif}})
stuff between endif and while

{{while a==1:}} infinite {{loop}} forever!
{{script}} a=0 {{endscript}}
{{for (a,b) in dictionary.items():}}
the key in the dictionary is {{a}} and the value is {{b}}.  And below is a script
{{script}}
        # THIS IS A SCRIPT
        x = 2
        y = 3
        # END OF THE SCRIPT
{{endscript}}
stuff after the script
{{endfor}}
stuff after the for stmt
{{endwhile}}
stuff after the while stmt

{{script}}
# test the free variables syntax error problem is gone
alpha = 3
def myfunction1(alpha=alpha):
    try:
        return free_variable # this would cause an error in an older version of preppy with python 2.2
    except:
        pass
    try:
        return alpha
    except:
        return "oops"
beta = myfunction1()
{{endscript}}
alpha = {{alpha}} and beta = {{beta}}

{${this is invalid but it's escaped, so no problem!}$}

end of text

{{script}}
# just a comment
{{endscript}}
stop here
"""

"""
# test code for quotestring
(qs, ds, c) = PreProcessor().quoteString(teststring, cursor=0)
print "---------quoted to ", c, `teststring[c:c+20]`
print qs
print "---------dict string"
print ds
"""
def dedent(text):
    """get rid of redundant indentation in text this dedenter IS NOT smart about converting tabs to spaces!!!"""
    lines = text.split("\n")
    # omit empty lines
    lempty = 0
    while lines:
        line0 = lines[0].strip()
        if line0 and line0[0]!='#': break
        del lines[0]
        lempty += 1
    if not lines: return "" # completely white
    line0 = lines[0]
    findfirstword = line0.find(line0.strip().split()[0])
    if findfirstword<0: raise ValueError('internal dedenting error')
    indent = line0[:findfirstword]
    linesout = []
    for l in lines:
        lines0 = l.strip()
        if not lines0 or lines0[0]=='#':
            linesout.append("")
            continue
        lindent = l[:findfirstword]
        if lindent!=indent:
            raise ValueError, "inconsistent indent expected %s got %s in %s" % (repr(indent), repr(lindent), l)
        linesout.append(l[findfirstword:])
    return '\n'.join(lempty*['']+linesout)


_pat = re.compile('{{\\s*|}}',re.M)
_s = re.compile(r'^(?P<start>while|if|elif|for)(?P<startend>\s+|$)|(?P<def>def\s*)(?P<defend>\(|$)|(?P<end>else|script|eval|endwhile|endif|endscript|endeval|endfor)(?:\s*$|(?P<endend>.+$))',re.DOTALL|re.M)

def _renumber(node,lineno_offset):
    if node.lineno: node.lineno += lineno_offset
    for child in node.getChildNodes(): _renumber(child,lineno_offset)

def _denumber(node,lineno=-1):
    if node.lineno!=lineno: node.lineno = lineno
    for child in node.getChildNodes(): _denumber(child,lineno)

class PreppyParser(pycodegen.Module):
    from compiler import parse as _cparse
    _cparse = staticmethod(_cparse)
    def __init__(self,*args,**kw):
        pycodegen.Module.__init__(self,*args,**kw)
        self._defSeen = 0

    def compile(self, display=0):
        tree = self._get_tree()
        gen = _preppy_ModuleCodeGenerator(tree)
        if display:
            import pprint
            print pprint.pprint(tree)
        self.code = gen.getCode()

    def getPycHeader(self):
        try:
            if getattr(self,'nosourcefile',0): raise ValueError
            mtime = os.path.getmtime(self.filename)
        except:
            import time
            mtime = time.time()
        mtime = struct.pack('<i', mtime)
        return self.MAGIC + mtime

    def _get_tree(self):
        from compiler import misc, syntax
        tree = self.__get_ast()
        misc.set_filename(self.filename, tree)
        syntax.check(tree)
        return tree

    def __lexerror(self, msg, pos):
        text = self.source
        pos0 = text.rfind('\n',0,pos)+1
        pos1 = text.find('\n',pos)
        if pos1<0: pos1 = len(text)
        raise SyntaxError('%s\n%s\n%s (near line %d of %s)' %(text[pos0:pos1],(' '*(pos-pos0)),msg,text.count('\n',0,pos)+1,self.filename))

    def __tokenize(self):
        text = self.source
        self._tokens = tokens = []
        a = tokens.append
        state = 0
        ix = 0
        for i in _pat.finditer(text):
            i0 = i.start()
            i1 = i.end()
            lineno = text.count('\n',0,ix)+1
            if i.group()!='}}':
                if state:
                    self.__lexerror('Unexpected {{', i0)
                else:
                    state = 1
                    if i0!=ix: a(('const',lineno,ix,i0))
                    ix = i1
            elif state:
                    state = 0
                    #here's where a preppy token is finalized
                    m = _s.match(text[ix:i0])
                    if m:
                        t = m.group('start')
                        if t:
                            if not m.group('startend'): self.__lexerror('Bad %s' % t, i0)
                            ix += len(t)    #skip over the 'while ' or 'for ' etc
                        else:
                            t = m.group('end')
                            if t:
                                ee = m.group('endend')
                                if ee and t!='else' and ee.strip()!=':': self.__lexerror('Bad %s' % t, i0)
                            else:
                                t = m.group('def')
                                if t:
                                    if not m.group('defend'): self.__lexerror('Bad %s' % t, i0)
                                    if self._defSeen:
                                        if self._defSeen>0:
                                            self.__lexerror('Only one def may be used',i0)
                                        else:
                                            self.__lexerror('def must come first',i0)
                                    else:
                                        self._defSeen = 1
                                ix += len(t)    #skip over the 'while ' or 'for ' etc
                    else:
                        t = 'expr'  #expression
                    if not self._defSeen: self._defSeen = -1
                    if i0!=ix: a((t,lineno,ix,i0))
                    ix = i1
        else:
            lineno = 0
        if state: self.__lexerror('Unterminated preppy token', ix)
        textLen = len(text)
        if ix!=textLen:
            lineno = text.count('\n',0,ix)+1
            a(('const',lineno,ix,textLen))
        a(('eof',lineno+1,textLen,textLen))
        self._tokenX = 0
        return tokens

    def __tokenText(self,colonRemove=0, strip=1):
        t = self._tokens[self._tokenX]
        text = self.source[t[2]:t[3]]
        if strip: text = text.strip()
        if colonRemove and text.endswith(':'): text = text[:-1]
        return unescape(text)

    def __tokenPop(self):
        t = self._tokens[self._tokenX]
        self._tokenX += 1
        return t

    def __preppy(self,funcs=['const','expr','while','if','for','script', 'eval','def'],followers=['eof']):
        from compiler.ast import Stmt
        C = []
        a = C.append
        mangle = '_%s__'%self.__class__.__name__
        tokens = self._tokens
        while 1:
            t = tokens[self._tokenX][0]
            if t in followers: break
            p = t in funcs and getattr(self,mangle+t) or self.__serror
            r = p()
            if type(r) is list: C += r
            elif r is not None: a(r)
        self.__tokenPop()
        return Stmt(C)
        
    def __def(self,followers=['endwhile']):
        from compiler.ast import While
        try:
            F = self._cparse('def X%s: pass' % self.__tokenText(colonRemove=1).strip()).node.nodes[0]
            self._fnc_argnames = F.argnames
            self._fnc_defaults = F.defaults
            self._fnc_varargs = F.varargs
            self._fnc_kwargs = F.kwargs
            self._fnc_flags = F.flags
        except:
            self.__error()
        t = self.__tokenPop()
        self._fnc_lineno = t[1]
        return None

    def __while(self,followers=['endwhile']):
        from compiler.ast import While
        try:
            cond = self._cparse(self.__tokenText(colonRemove=1),'eval').node
        except:
            self.__error()
        t = self.__tokenPop()
        _renumber(cond,t[1]-1)
        r = While(cond,self.__preppy(followers=followers),None)
        return r

    def __for(self,followers=['endfor']):
        from compiler.ast import For
        try:
            f = self._cparse('for %s:\n pass\n'%self.__tokenText(colonRemove=1),'exec').node.nodes[0]
        except:
            self.__error()
        t = self.__tokenPop()
        _renumber(f,t[1]-1)
        f.body = self.__preppy(followers=followers)
        return f

    def __script(self,mode='script'):
        self.__tokenPop()
        text = dedent(self.__tokenText(strip=0))
        scriptMode = 'script'==mode
        if text:
            try:
                stmt = self._cparse(text,scriptMode and 'exec' or 'eval').node
            except:
                self.__error()
        t = self.__tokenPop()
        end = 'end'+mode
        try:
            assert self._tokens[self._tokenX][0]==end
            self.__tokenPop()
        except:
            self.__error(end+' expected')
        if not text: return []
        _renumber(stmt,t[1]-1)

        if scriptMode: return stmt.nodes
        from compiler.ast import Discard, CallFunc, Name, Const
        return Discard(CallFunc(Name('__swrite__'), [stmt.getChildren()[0]],None,None))

    def __eval(self):
        return self.__script(mode='eval')

    def __if(self,followers=['endif','elif','else']):
        from compiler.ast import If
        tokens = self._tokens
        CS = []
        t = 'elif'
        while t=='elif':
            try:
                cond = self._cparse(self.__tokenText(colonRemove=1),'eval').node
            except:
                self.__error()
            t = self.__tokenPop()
            _renumber(cond,t[1]-1)
            CS.append((cond,self.__preppy(followers=followers)))
            t = tokens[self._tokenX-1][0]   #we consumed the terminal in __preppy
            if t=='elif': self._tokenX -= 1
        if t=='else':
            stmt = self.__preppy(followers=['endif'])
        else:
            stmt = None
        return If(CS,stmt)

    def __const(self):
        from compiler.ast import Discard, CallFunc, Name, Const
        try:
            n = Discard(CallFunc(Name('__write__'), [Const(self.__tokenText(strip=0))], None, None))
            t = self.__tokenPop()
            _renumber(n,t[1]-1)
            return n
        except:
            self.__error('bad constant')

    def __expr(self):
        from compiler.ast import Discard, CallFunc, Name
        t = self.__tokenText()
        try:
            n= Discard(CallFunc(Name('__swrite__'), [self._cparse(t,'eval').getChildren()[0]],None,None))
            t = self.__tokenPop()
            _renumber(n,t[1]-1)
            return n
        except:
            self.__error('bad expression')

    def __error(self,msg='invalid syntax'):
        pos = self._tokens[self._tokenX][2]
        import traceback, StringIO
        f = StringIO.StringIO()
        traceback.print_exc(file=f)
        f = f.getvalue()
        m = 'File "<string>", line '
        i = f.rfind('File "<string>", line ')
        if i<0:
            t, v = map(str,sys.exc_info()[:2])
            self.__lexerror('%s %s(%s)' % (msg, t, v),pos)
        else:
            i += len(m)
            f = f[i:].split('\n')
            n = int(f[0].strip())+self.source[:pos].count('\n')
            raise SyntaxError('  File %s, line %d\n%s' % (self.filename,n,'\n'.join(f[1:])))

    def __serror(self,msg='invalid syntax'):
        self.__lexerror(msg,self._tokens[self._tokenX][2])

    def __parse(self,text=None):
        from compiler.ast import Stmt
        if text: self.source = text
        self.__tokenize()
        return self.__preppy().nodes

    def __get_ast(self):
        from compiler.ast import    Module, Stmt, Assign, AssName, Const, Function, For, Getattr,\
                                    TryFinally, TryExcept, If, Import, AssAttr, Name, CallFunc,\
                                    Class, Compare, Raise, And, Mod, Tuple, Pass, Not, Exec, List,\
                                    Discard, Keyword, Return, Dict, Break, AssTuple, Subscript,\
                                    Printnl, From, Lambda

        preppyNodes = self.__parse()
        if self._defSeen==1:
            fixargs = self._fnc_argnames
            defaults = list(self._fnc_defaults)
            if self._fnc_kwargs:
                spargs = [fixargs[-1]]
                fixargs = fixargs[:-1]
            else:
                spargs = ['__kwds__']
            if self._fnc_varargs:
                spargs.insert(0,fixargs[-1])
                fixargs = fixargs[:-1]
            kwargs = fixargs[-len(defaults):]
            fixargs = fixargs[:-len(defaults)]
            flags = self._fnc_flags

            #construct the getOutput function
            nodes = [Assign([AssName('__lquoteFunc__', 'OP_ASSIGN')], CallFunc(Getattr(Name(spargs[-1]), 'setdefault'), [Const('__lquoteFunc__'), Name('str')], None, None)),
                    Discard(CallFunc(Getattr(Name(spargs[-1]), 'pop'),[Const('__lquoteFunc__')], None, None)),
                    Assign([AssName('__quoteFunc__', 'OP_ASSIGN')], CallFunc(Getattr(Name(spargs[-1]), 'setdefault'), [Const('__quoteFunc__'), Name('str')], None, None)),
                    Discard(CallFunc(Getattr(Name(spargs[-1]), 'pop'), [Const('__quoteFunc__')], None, None))]
            if not self._fnc_kwargs:
                nodes += [If([(Name(spargs[-1]), Stmt([Raise(CallFunc(Name('TypeError'), [Const('get: unexpected keyword arguments')], None, None), None, None)]))], None)]
            nodes += [Assign([AssName('__append__', 'OP_ASSIGN')], Getattr(List(()), 'append')),
                    Assign([AssName('__write__', 'OP_ASSIGN')], Lambda(['x'], [], 0, CallFunc(Name('__append__'), [CallFunc(Name('__lquoteFunc__'), [Name('x')], None, None)], None, None))),
                    Assign([AssName('__swrite__', 'OP_ASSIGN')], Lambda(['x'], [], 0, CallFunc(Name('__append__'), [CallFunc(Name('__quoteFunc__'), [Name('x')], None, None)], None, None)))]
            for n in nodes: _denumber(n,self._fnc_lineno)
            preppyNodes = nodes + preppyNodes + [Return(CallFunc(Getattr(Const(''), 'join'), [Getattr(Name('__append__'), '__self__')], None, None))]
            argnames = list(fixargs)+list(kwargs)+list(spargs)
            FA = ('get',argnames, defaults,flags|8,None,Stmt(preppyNodes))
            global _newPreambleAst
            if not _newPreambleAst:
                _newPreambleAst = self._cparse(_newPreamble).node.nodes
                map(_denumber,_newPreambleAst)
            extraAst = _newPreambleAst
        else:
            global _preambleAst, _localizer
            if not _preambleAst:
                _preambleAst = self._cparse(_preamble).node.nodes
                map(_denumber,_preambleAst)
                _localizer = [Assign([AssName('__d__', 'OP_ASSIGN')], Name('dictionary')), Discard(CallFunc(Getattr(CallFunc(Name('globals'), [], None, None), 'update'), [Name('__d__')], None, None))]
            preppyNodes = _localizer+preppyNodes
            FA = ('__code__', ['dictionary', 'outputfile', '__write__','__swrite__','__save_sys_stdout__'], (), 0, None, Stmt(preppyNodes))
            extraAst = _preambleAst
        if sys.hexversion >=0x2040000: FA = (None,)+FA
        return Module(self.filename,
                Stmt([Assign([AssName('__checksum__', 'OP_ASSIGN')], Const(getattr(self,'sourcechecksum'))),
                        Function(*FA),
                    ]+extraAst))

_preambleAst=None
_preamble='''def run(dictionary, __write__=None, quoteFunc=str, outputfile=None,code=__code__):
    ### begin standard prologue
    import sys
    __save_sys_stdout__ = sys.stdout
    try: # compiled logic below
        try:
            # if outputfile is defined, blindly assume it supports file protocol, reset sys.stdout
            if outputfile is None:
                raise NameError
            stdout = sys.stdout = outputfile
        except NameError:
            stdout = sys.stdout
            outputfile = None
        # make sure quoteFunc is defined:
        if quoteFunc is None:
            raise ValueError("quoteFunc must be defined")
        # make sure __write__ is defined
        try:
            if __write__ is None:
                raise NameError
            if outputfile and __write__:
                raise ValueError, "do not define both outputfile (%s) and __write__ (%s)." %(outputfile, __write__)
            class stdout: pass
            stdout = sys.stdout = stdout()
            stdout.write = lambda x:__write__(quoteFunc(x))
        except NameError:
            __write__ = lambda x: stdout.write(quoteFunc(x))
        __swrite__ = lambda x: __write__(quoteFunc(x))
        code(dictionary,outputfile,__write__,__swrite__,__save_sys_stdout__)
    finally: #### end of compiled logic, standard cleanup
        import sys # for safety
        #print "resetting stdout", sys.stdout, "to", __save_sys_stdout__
        sys.stdout = __save_sys_stdout__
def getOutput(dictionary, quoteFunc=str):
    buf=[]
    run(dictionary,__write__=buf.append, quoteFunc=quoteFunc)
    return ''.join(buf)

def getOutputFromKeywords(quoteFunc=str, **kwds):
    buf=[]
    run(kwds,__write__=buf.append, quoteFunc=quoteFunc)
    return ''.join(buf)

if __name__=='__main__':
    run()
'''
_newPreambleAst=None
_newPreamble='''def run(*args,**kwds):
    raise ValueError('Wrong kind of prep file')
def getOutput(*args,**kwds):
    run()
if __name__=='__main__':
    run()
'''



def testgetOutput(name="testoutput"):
    mod = getModule(name,'.',savePyc=1,sourcetext=teststring,importModule=1)
    print mod.getOutput({})

def testgetmodule(name="testoutput"):
    #name = "testpreppy"
    print "trying to load", name
    result = getPreppyModule(name, verbose=1)
    print "load successful! running result"
    print "=" * 100
    result.run({})
try:
    from reportlab.lib.utils import rl_get_module
except:
    def rl_get_module(name,dir):
        f, p, desc= imp.find_module(name,[dir])
        if sys.modules.has_key(name):
            om = sys.modules[name]
            del sys.modules[name]
        else:
            om = None
        try:
            try:
                return imp.load_module(name,f,p,desc)
            except:
                raise ImportError('%s[%s]' % (name,dir))
        finally:
            if om: sys.modules[name] = om
            del om
            if f: f.close()

# cache found modules by source file name
FILE_MODULES = {}
SOURCE_MODULES = {}
def getModule(name,
              directory=".",
              source_extension=".prep",
              verbose=0,
              savefile=None,
              sourcetext=None,
              savePy=0,
              force=0,
              savePyc=1,
              importModule=1,_globals=None):
    """Returns a python module implementing the template, compiling if needed.

    force: ignore up-to-date checks and always recompile.
    """
    if isinstance(name,unicode): name = name.encode('utf8')
    if isinstance(directory,unicode): directory = directory.encode('utf8')
    if isinstance(source_extension,unicode): source_extension = source_extension.encode('utf8')
    if isinstance(sourcetext,unicode): sourcetext = sourcetext.encode('utf8')
    if hasattr(name,'read'):
        sourcetext = name.read()
        name = getattr(name,'name',None)
        if not name: name = '_preppy_'+md5(sourcetext+repr(VERSION)).digest()
    else:
        # it's possible that someone could ask for
        #  name "subdir/spam.prep" in directory "/mydir", instead of
        #  "spam.prep" in directory "/mydir/subdir".  Failing to get
        # this right means getModule can fail to find it and recompile
        # every time.  This is common during batch compilation of directories.
        extraDir, name = os.path.split(name)
        if extraDir:
            if os.path.isabs(extraDir):
                directory = extraDir
            else:
                directory = directory + os.sep + extraDir
        dir = os.path.abspath(os.path.normpath(directory))

        # they may ask for 'spam.prep' instead of just 'spam'.  Trim off
        # any extension
        name = os.path.splitext(name)[0]
        if verbose:
            print 'checking %s...' % os.path.join(dir, name),
        # savefile is deprecated but kept for safety.  savePy and savePyc are more
        # explicit and are the preferred.  By default it generates a pyc and no .py
        # file to reduce clutter.
        if savefile and savePyc == 0:
            savePyc = 1

    if sourcetext is not None:
        # they fed us the source explicitly
        if not name: name = '_preppy_'+md5(sourcetext+repr(VERSION)).digest()
        if verbose: print "sourcetext provided...",
        sourcefilename = "<input text %s>" % name
        sourcechecksum = md5(sourcetext + repr(VERSION)).digest()
        nosourcefile = 1
        module = SOURCE_MODULES.get(sourcetext,None)
        if module:
            return module
    else:
        nosourcefile = 0
        # see if the module exists as a python file
        sourcefilename = os.path.join(dir, name+source_extension)
        module = FILE_MODULES.get(sourcefilename,None)
        if module:
            return module

        try:
            module = rl_get_module(name,dir)
            module.__dict__['include'] = include
            if _globals:
                module.__dict__.update(_globals)
            checksum = module.__checksum__
            if verbose: print "found...",
        except: # ImportError:  #catch ALL Errors importing the module (eg name="")
            module = checksum = None
            if verbose: print " py/pyc not found...",
            # check against source file
        try:
            sourcefile = open(sourcefilename, "r")
        except:
            if verbose: print "no source file, reuse...",
            if module is None:
                raise ValueError, "couldn't find source %s or module %s" % (sourcefilename, name)
            # use the existing module??? (NO SOURCE PRESENT)
            FILE_MODULES[sourcefilename] = module
            return module
        else:
            sourcetext = sourcefile.read()
            # NOTE: force recompile on each new version of this module.
            sourcechecksum = md5(sourcetext + repr(VERSION)).digest()
            if sourcechecksum==checksum:
                if force==0:
                    # use the existing module. it matches
                    if verbose:
                        print "up to date."
                    FILE_MODULES[sourcefilename] = module
                    return module
                else:
                    # always recompile
                    if verbose: print 'forced recompile,',
            elif verbose:
                print "changed,",

    # if we got here we need to rebuild the module from source
    if verbose: print "recompiling"
    global DIAGNOSTIC_FUNCTION
    DIAGNOSTIC_FUNCTION = None
    P = PreppyParser(sourcetext,sourcefilename)
    P.nosourcefile = nosourcefile
    P.sourcechecksum = sourcechecksum
    P.compile(0)

    # default is compile to bytecode and save that.
    if savePyc:
        f = open(dir + os.sep + name + '.pyc','wb')
        P.dump(f)
        f.close()

    # now make a module
    from imp import new_module
    module = new_module(name)
    module.__dict__['include'] = include
    if _globals:
        module.__dict__.update(_globals)
    exec P.code in module.__dict__
    if importModule:
        if nosourcefile:
            SOURCE_MODULES[sourcetext] = module
        else:
            FILE_MODULES[sourcefilename] = module
    return module

# support the old form here
getPreppyModule = getModule

####################################################################
#
#   utilities for compilation, setup scripts, housekeeping etc.
#
####################################################################
def installImporter():
    "This lets us import prep files directly"
    # the python import mechanics are only invoked if you call this,
    # since preppy has very few dependencies and I don't want to
    #add to them.
    from ihooks import ModuleLoader, ModuleImporter, install
    class PreppyLoader(ModuleLoader):
        "This allows prep files to be imported."

        def find_module_in_dir(self, name, dir, allow_packages=1):
            ModuleLoader = self.__class__.__bases__[0]
            stuff = ModuleLoader.find_module_in_dir(self, name, dir, allow_packages)
            if stuff:
                #print 'standard module loader worked'
                return stuff
            else:
                if dir:
                    prepFileName = dir + os.sep + name + '.prep'
                else:
                    prepFileName = name + '.prep'

                if os.path.isfile(prepFileName):
                    #compile WITHOUT IMPORTING to avoid triggering recursion
                    mod = compileModule(prepFileName, verbose=0, importModule=0)
                    #now use the default...
                    return ModuleLoader.find_module_in_dir(self, name, dir, allow_packages)
                else:
                    return None

    loader = PreppyLoader()
    importer = ModuleImporter(loader=loader)
    install(importer)

def uninstallImporter():
    import ihooks
    ihooks.uninstall()

def compileModule(fn, savePy=0, force=0, verbose=1, importModule=1):
    "Compile a prep file to a pyc file.  Optionally, keep the python source too."
    name, ext = os.path.splitext(fn)
    d = os.path.dirname(fn)
    return getModule(os.path.basename(name), directory=d, source_extension=ext,
                     savePyc=1, savePy=savePy, force=force,
                     verbose=verbose, importModule=importModule)

def compileModules(pattern, savePy=0, force=0, verbose=1):
    "Compile all prep files matching the pattern."
    import glob
    filenames = glob.glob(pattern)
    for filename in filenames:
        compileModule(filename, savePy, force, verbose)

from fnmatch import fnmatch
def compileDir(dirName, pattern="*.prep", recursive=1, savePy=0, force=0, verbose=1):
    "Compile all prep files in directory, recursively if asked"
    if verbose: print 'compiling directory %s' % dirName
    if recursive:
        def _visit(A,D,N,pattern=pattern,savePy=savePy, verbose=verbose,force=force):
            for filename in filter(lambda fn,pattern=pattern: fnmatch(fn,pattern),
                    filter(os.path.isfile,map(lambda n, D=D: os.path.join(D,n),N))):
                compileModule(filename, savePy, force, verbose)
        os.path.walk(dirName,_visit,None)
    else:
        compileModules(os.path.join(dirName, pattern), savePy, force, verbose)

def _cleanFiles(filenames,verbose):
    for filename in filenames:
        if verbose:
            print '  found ' + filename + '; ',
        root, ext = os.path.splitext(os.path.abspath(filename))
        done = 0
        if os.path.isfile(root + '.py'):
            os.remove(root + '.py')
            done = done + 1
            if verbose: print ' removed .py ',
        if os.path.isfile(root + '.pyc'):
            os.remove(root + '.pyc')
            done = done + 1
            if verbose: print ' removed .pyc ',
        if done == 0:
            if verbose:
                print 'nothing to remove',
        print

def cleanDir(dirName, pattern="*.prep", recursive=1, verbose=1):
    "Removes all py and pyc files matching any prep files found"
    if verbose: print 'cleaning directory %s' % dirName
    if recursive:
        def _visit(A,D,N,pattern=pattern,verbose=verbose):
            _cleanFiles(filter(lambda fn,pattern=pattern: fnmatch(fn,pattern),
                    filter(os.path.isfile,map(lambda n, D=D: os.path.join(D,n),N))),verbose)
        os.path.walk(dirName,_visit,None)
    else:
        import glob
        _cleanFiles(filter(os.path.isfile,glob.glob(os.path.join(dirName, pattern))),verbose)

def compileStuff(stuff, savePy=0, force=0, verbose=0):
    "Figures out what needs compiling"
    #print 'compileStuff...savePy=%d, force=%d, verbose=%d' % (savePy, force, verbose)
    if os.path.isfile(stuff):
        compileModule(stuff, savePy=savePy, force=force, verbose=verbose)
    elif os.path.isdir(stuff):
        compileDir(stuff, savePy=savePy, force=force, verbose=verbose)
    else:
        compileModules(stuff, savePy=savePy, force=force, verbose=verbose)

def extractKeywords(arglist):
    "extracts a dictionary of keywords"
    d = {}
    for arg in arglist:
        chunks = arg.split('=')
        if len(chunks)==2:
            key, value = chunks
            d[key] = value
    return d

from reportlab.lib.utils import find_locals
def _find_quoteValue(name):
    def func(L):
        if L.has_key(name): return L[name] 
        return None
    try:
        return find_locals(func)
    except:
        return None

def include(viewName,*args,**kwd):
    dir, filename = os.path.split(viewName)
    root, ext = os.path.splitext(filename)
    m = {}
    if dir: m['directory'] = dir
    if ext: m['source_extension'] = ext
    m = getModule(root,**m)
    if hasattr(m,'get'):
        #newstyle
        lquoter = quoter = None
        if kwd.has_key('__quoteFunc__'):
            quoter = kwd.pop('__quoteFunc__')
        elif kwd.has_key('quoteFunc'):
            quoter = kwd.pop('quoteFunc')
        if not quoter:
            quoter = _find_quoteValue('__quoteFunc__')
            if not quoter:
                quoter = _find_quoteValue('quoteFunc')
        if kwd.has_key('__lquoteFunc__'):
            lquoter = kwd.pop('__lquoteFunc__')
        if not lquoter:
            lquoter = _find_quoteValue('__lquoteFunc__')
            if not lquoter:
                lquoter = _find_quoteValue('quoteFunc')
        return m.get(__quoteFunc__=quoter or str,__lquoteFunc__=lquoter or str, *args,**kwd)
    else:
        #oldstyle
        if args:
            if len(args)>1:
                raise TypeError("include for old style prep file can have only one positional argument, dictionary")
            if kwd.has_key('dictionary'):
                raise TypeError('include: dictionary argument specified twice')
            dictionary = args(1).copy()
        elif kwd.has_key('dictionary'):
            dictionary = kwd.pop('dictionary').copy()
        else:
            dictionary = {}
        quoteFunc = None
        if kwd.has_key('quoteFunc'):
            quoteFunc = kwd.pop('quoteFunc')
        elif kwd.has_key('__quoteFunc__'):
            quoteFunc = kwd.pop('__quoteFunc__')
        if not quoteFunc:
            quoteFunc = _find_quoteValue('quoteFunc')
            if not quoteFunc:
                quoteFunc = _find_quoteValue('__quoteFunc__')
        dictionary.update(kwd)
        return m.getOutput(dictionary,quoteFunc=quoteFunc)

def main():
    if len(sys.argv)>1:
        name = sys.argv[1]

        if name == 'compile':
            names = sys.argv[2:]

            # save the intermediate python file
            if '--savepy' in names:
                names.remove('--savepy')
                savePy = 1
            elif '-p' in names:
                names.remove('-p')
                savePy = 1
            else:
                savePy = 0

            # force recompile every time
            if '--force' in names:
                names.remove('--force')
                force = 1
            elif '-f' in names:
                names.remove('-f')
                force = 1
            else:
                force = 0

            # extra output
            if '--verbose' in names:
                names.remove('--verbose')
                verbose = 1
            elif '-v' in names:
                names.remove('-v')
                verbose = 1
            else:
                verbose = 0

            for arg in names:
                compileStuff(arg, savePy=savePy, force=force, verbose=verbose)

        elif name == 'clean':
            for arg in sys.argv[2:]:
                cleanDir(arg, verbose=1)

        elif name == 'run':
            moduleName = sys.argv[2]
            params = extractKeywords(sys.argv)
            module = getPreppyModule(moduleName, verbose=0)
            module.run(params)

        else:
            #default is run
            moduleName = sys.argv[1]
            module = getPreppyModule(moduleName, verbose=0)
            if hasattr(module,'get'):
                print module.get()
            else:
                params = extractKeywords(sys.argv)
                module.run(params)
    else:
        print "no argument: running tests"
        testgetOutput()
        print; print "PAUSING.  To continue hit return"
        raw_input("now: ")
        testgetmodule()

if __name__=="__main__":
    main()
