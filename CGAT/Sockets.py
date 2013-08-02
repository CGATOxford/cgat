################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
"""Sockets.py - working with sockets
====================================

This class allows you to send variable length strings
over a socket. As I am new at this, it is probably not
efficient and 100% fault tolerant.

Look at the end of the file for example usage.
"""

import sys, os, string, getopt

import socket

class SocketException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
    
class Socket:

    mBufferSize = 4096
    mHeaderSize = 10

    def __init__( self, sock = None):

        if sock:
            self.mSocket = sock
        else:
            self.mSocket = socket.socket( socket.AF_INET, socket.SOCK_STREAM)

    def __del__(self):
        self.Close()

    def MakeServer( self, host, port ):
        hostname = socket.gethostbyname( host )
        self.mSocket.bind( (hostname, port))
        self.mSocket.listen(1)

    def MakeClient( self, host, port ):
        hostname = socket.gethostbyname( host )
        self.mSocket = socket.socket( socket.AF_INET, socket.SOCK_STREAM)
        self.mSocket.connect( (hostname, port) )

    def SendMessage( self, msg ):
        """send atribraty length message over socket."""
        
        lmsg = len(msg)
    
        m = str(lmsg) + " " * (self.mHeaderSize - len(str(lmsg)))
        
        sent = self.mSocket.send( m )
        
        if sent != self.mHeaderSize: raise SocketException("msg header could not be sent")
        
        sent = 0
        while sent < lmsg:
            r = self.mSocket.send( msg[sent:sent+self.mBufferSize] )
            if r == 0: raise SocketException("connection broken")
            sent += r
            
    def ReceiveMessage( self ):
        """receive arbitrary length message.
        """

        msg = ""
        while len(msg) < self.mHeaderSize:
            chunk = self.mSocket.recv( self.mHeaderSize - len(msg)  )
            if chunk == "": raise SocketException("connection broken")
            msg += chunk
        
        lmsg = string.atoi(msg)
        msg = ""
        
        while lmsg > 0:
            chunk = self.mSocket.recv( min(self.mBufferSize, lmsg) )
            if chunk == "": raise SocketException("connection broken")
            msg += chunk
            lmsg -= len(chunk)
            
        return msg

    def Accept( self ):
        """wait for connection to come."""
        conn, addr = self.mSocket.accept()
        return Socket( conn ), addr

    def Close( self ):
        """close socket."""
        self.mSocket.close()
        
if __name__ == "__main__":

    param_mode = None
    param_host = "kerberos.biocenter.helsinki.fi"
    param_port = 9000
    param_loglevel = 4

    try:
        optlist, args = getopt.getopt(sys.argv[1:],
                                      "h:",
                                      ["host="])
        
    except getopt.error, msg:
        print USAGE
        print msg
        sys.exit(2)

    for o,a in optlist:
        if o in ("-h", "--host"):
            param_host = a
    
    if len(args) != 1:
        raise "not enough arguments specified"

    param_mode = args[0]

    if param_loglevel >= 2:
        print "# host=", param_host
        print "# port=", param_port
        print "# mode=", param_mode
        
    sock = Socket()
    
    if param_mode == "server":
        
        print "gtg_align server starting"
        sock.MakeServer( param_host, param_port )

        while 1:
            
            conn, addr = sock.Accept()
            
            msg = conn.ReceiveMessage()
            print "server received", msg

            msg += "a"
            
            conn.SendMessage( msg)
            print "server sent", msg
            
        sock.Close()

    elif param_mode == "client":
        print "client starting"        
        msg = "b"

        for x in range(1, 10):

            sock.MakeClient( param_host, param_port )

            sock.SendMessage( msg)
            print "client sent", msg

            msg = sock.ReceiveMessage()
            print "client received", msg
            
            msg += "b"
            
            sock.Close()

            
