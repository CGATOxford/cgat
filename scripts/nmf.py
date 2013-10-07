'''
nmf.py - compute NMF factorization
==================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Needs to adapted as a libary.

Usage
-----

Example::

   python nmf.py --help

Type::

   python nmf.py --help

for command line help.

Command line options
--------------------

'''
from numpy import *
from numpy.linalg import norm
from time import time
from sys import stdout

def nmf(V,Winit,Hinit,tol,timelimit,maxiter,log = None):
   """
   (c,W,H) = nmf(V,Winit,Hinit,tol,timelimit,maxiter)
   W,H: output solution
   c: output flag to indicate convergence
   Winit,Hinit: initial solution
   tol: tolerance for a relative stopping condition
   timelimit, maxiter: limit of time and iterations
   log: logfile
   """

   W = Winit
   H = Hinit
   initt = time()

   gradW = dot(W, dot(H, H.T)) - dot(V, H.T)
   gradH = dot(dot(W.T, W), H) - dot(W.T, V)
   initgrad = norm(r_[gradW, gradH.T])

   if log: log.write('# init gradient norm %f\n# progress: ' % initgrad )

   tolW = max(0.001,tol)*initgrad
   tolH = tolW

   converged = False
   for iter in xrange(1,maxiter):
      # stopping condition
      projnorm = norm(r_[gradW[logical_or(gradW<0, W>0)],
                         gradH[logical_or(gradH<0, H>0)]])

      if projnorm < tol*initgrad: 
         if log: log.write( "\n# converged\n" )
         converged = True
         break

      if time() - initt > timelimit: 
         if log: log.write( "\n# time limit reached\n" )
         converged = False
         break

      (W, gradW, iterW) = nlssubprob(V.T,H.T,W.T,tolW,1000, log)
      W = W.T
      gradW = gradW.T

      if iterW==1: tolW = 0.1 * tolW

      (H,gradH,iterH) = nlssubprob(V,W,H,tolH,1000, log)
      if iterH==1: tolH = 0.1 * tolH

      if log and iter % 10 == 0: 
         log.write('.')
         log.flush()
   else:
      converged = False
      if log: log.write( "\n# iteration limit reached\n" )

   if log: log.write( '# Iter = %d Final proj-grad norm %f\n' % (iter, projnorm) )
   return (converged,W,H)

def nlssubprob(V,W,Hinit,tol,maxiter, log = None):
   """
   H, grad: output solution and gradient
   iter: #iterations used
   V, W: constant matrices
   Hinit: initial solution
   tol: stopping tolerance
   maxiter: limit of iterations
   """
 
   H = Hinit
   WtV = dot(W.T, V)
   WtW = dot(W.T, W) 

   alpha = 1; beta = 0.1;
   for iter in xrange(1, maxiter):  
      grad = dot(WtW, H) - WtV
      projgrad = norm(grad[logical_or(grad < 0, H >0)])
      if projgrad < tol: break

      # search step size 
      for inner_iter in xrange(1,20):
         Hn = H - alpha*grad
         Hn = where(Hn > 0, Hn, 0)
         d = Hn-H
         gradd = sum(grad * d)
         dQd = sum(dot(WtW,d) * d)
         suff_decr = 0.99*gradd + 0.5*dQd < 0
         if inner_iter == 1:
            decr_alpha = not suff_decr
            Hp = H
         if decr_alpha: 
            if suff_decr:
               H = Hn
               break
            else:
               alpha = alpha * beta;
         else:
            if not suff_decr or (Hp == Hn).all():
               H = Hp
               break
            else:
               alpha = alpha/beta
               Hp = Hn

      if iter == maxiter:
         if log: log.write( '# max iteration in nlssubprob\n' )
   return (H, grad, iter)
