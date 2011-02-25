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
Database.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import os, sys, time, string, traceback, re
import MySQLdb
import _mysql           # necessary for exception handling

WARNINGS = 1                    # print out warnings?

class Database:
    dbhandle = None

    socket = None
    host = None
    user = None
    passwd = None
    port = None
    dbname = None

    ## perform retries if you lost connection
    ## by getting a new connection.
    mRetryLostConnection = 1

    
    mRetryTimeIntervall  = 60


    def __init__(self):
        pass

    def AddOptions( self, parser ):
        """add connection options to a command line parser.
        """
        
        parser.add_option("-D", "--Database", dest="database", type="string",
                          help="database to connect to.")
        
        parser.add_option("-H", "--Host", dest="host", type="string",
                          help="host to connect to.")

        parser.add_option("-U", "--Username", dest="user", type="string",
                          help="username in database.")

        parser.add_option("-W", "--Password", dest="passwd", type="string",
                          help="password for database.")

        parser.add_option("-S", "--Socket", dest="socket", type="string",
                          help="socket of database connection.")

        parser.add_option("-P", "--Port", dest="port", type="string",
                          help="port to use.")

        parser.set_defaults(
            host = None,
            database = "",
            user = None,
            passwd = "",
            port = None, 
            socket = "/tmp/mysql.sock" )
            
    def Connect( self, options ):

	try:
            if options.port:
                self.dbhandle = MySQLdb.connect(host        = options.host,
                                                user        = options.user,
                                                passwd      = options.password,
                                                db          = options.database,
                                                port        = options.port)
            else:
                self.dbhandle = MySQLdb.connect(host        = options.host,
                                                user        = options.user,
                                                passwd      = options.password,
                                                db          = options.database,
                                                unix_socket = options.socket )
	    return self.dbhandle
	except _mysql.OperationalError, detail:
            print "------------------------"
	    print "MySQL raised error %s: %s!" % detail.args
            traceback.print_stack()
            print "------------------------"
	    return 0

    def Close( self ):
        """close connection. No further activity possible.
        """
        self.dbhandle.close()
	    
    def Execute( self, statement ):
        """execute a mysql statement and return a cursor object.

        Die, if there is an error.
        """
	try:
	    c = self.dbhandle.cursor()
	    c.execute( statement )
        except _mysql.ProgrammingError, detail:
            
            print """# mysql raised programming error %s: %s!
# %s""" % (detail.args[0], detail.args[1], statement)
            traceback.print_stack(file=sys.stdout)                        
            sys.exit(1)
            
	except _mysql.OperationalError, detail:
            print """# mysql raised operational error %s: %s!
# %s""" % (detail.args[0], detail.args[1], statement)
            traceback.print_stack(file=sys.stdout)                        
            
            if detail.args[0] == "2013" and self.mRetryLostConnection:
            
                while 1:
            
                    sys.stderr.write( "retrying to establish connection\n" )
                    sys.stderr.flush()
                
                    time.sleep(self.mRetryTimeIntervall)
                    if self.Connect( host = self.host,
                                     user = self.user,
                                     passwd = self.passwd,
                                     dbname = self.dbname,
                                     port = self.port,
                                     socket = self.socket ):
                        break

                ## carefull: infinite loop!
                return self.Execute(statement)
            
            sys.exit(1)
                
	except _mysql.Warning, detail:
            if WARNINGS:
                print "------------------------"
                print "MySQL raised a warning: %s!" % detail.args
                traceback.print_stack(file=sys.stdout)
                print "Offending statement:"
                print statement
                print "------------------------"

        ## check, whether query was interrupted:
        errno = self.dbhandle.errno()
        if errno != 0:
            print """mysql error %i: %s!
%s""" % ( errno, self.dbhandle.error(), statement)
            traceback.print_stack(file=sys.stdout)
            sys.exit(1)
            
        return c

    def QuoteString( self, text ):
	return _mysql.escape_string(text)

    def GetDate(self):			# return date in format yyyy-mm-dd
	return time.strftime("%Y-%m-%d", time.localtime(time.time()))

    def Create( self ):
        if len(self.socket) == 0:
            self.socket = "/tmp/mysql.sock"             # use standard socket under unix

	try:
	    self.dbhandle = MySQLdb.connect(host        = self.host,
					    user        = self.user,
					    passwd      = self.passwd,
					    unix_socket = self.socket )
	except _mysql.OperationalError, detail:
	    print "MySQL raised error %s: %s!" % detail.args  
	    return 0

	c  = self.dbhandle.cursor()
	
	try:
	    c.execute('CREATE DATABASE %s' % self.dbname)
            return 1
	except _mysql.OperationalError, detail:
	    print "MySQL raised error %s: %s!" % detail.args
            return 0
    
    def Get_Tables( self, database = None):
        
        if database:
            statement = "SHOW TABLES FROM " + database
        else:
            statement = "SHOW TABLES"
            
        return map(lambda x: x[0], self.Execute(statement).fetchall())

    def DropTable( self, table_name):
        """drop table.
        """
        return self.Execute("DROP TABLE IF EXISTS %s" % table_name)

    def Exists( self, table_name ):
        """returns 1, if a table with name name exists."""
        # How can I locate an element in a list without iterating through it
        # myself? The code below looks ugly (exceptions, but this is the Python style?).

        r = string.split(table_name, ".")
        if len(r) > 1:
            database = r[0]
            table = r[1]
        else:
            database = None
            table = r[0]
        
        tables = self.Get_Tables(database)
        try:
            tables.index(table)
            return 1
        except (ValueError):
            return 0


    #--------------------------------------------------------------------------------------
    def GetDatabase( self ):
        """return current database name."""
        return self.dbname
    
    #--------------------------------------------------------------------------------------    
    def UseDatabase( self, db_name ):
        """switch database."""
        self.dbname = db_name
        return self.Execute( "USE %s " % db_name )

    #--------------------------------------------------------------------------------------
    def GetUser(self):
        return self.user
    #--------------------------------------------------------------------------------------
    def GetHost(self):
        return self.host
    #--------------------------------------------------------------------------------------
    def GetPort(self):
        return self.port
    #--------------------------------------------------------------------------------------
    def GetPassword(self):
        return self.passwd
    
    #------------------------------------------------------------------------------------------------------
    def SelectIntoOutfile( self, statement, filename ):
        """Redirect SELECT-statement INTO outfile on local disc. As this
        is not directly possible via mysqld, mysql -e > outfile is used.

        Note: password is not given here, it has to be set in .my.cnf,
        as otherwise it is visible via ps to all.

        options for mysql:
        -N: skip column headers
        """
        
        if os.system( "mysql -N -h%s -u%s -P%i -e '%s' > %s" % (self.host, self.user,
                                                                self.port,
                                                                statement, filename)):
            print "--> error causing statement: %s" % statement
        
    #------------------------------------------------------------------------------------------------------
    def LoadFile( self, filename, statement ):
        """load a file using statement.
        """
        
        os.chmod( filename, 0664)
        self.Execute( statement % filename )
        
        return 

def executewait( dbhandle, statement, error = None, retries = -1, wait=5):
    '''execute an sql ``statement``.

    If ``error`` is given, it is scanned for locking issues.

    Retry ``retries`` times if set to a positive number.
    A retry of ``0`` indicates no retry, a negative number retries
    infinitely. 

    The process waits ``wait`` seconds between each retry.

    Returns a cursor object.
    '''
    if error == None: error = Exception

    cc = dbhandle.cursor()    

    while 1:
        try:
            cc.execute( statement )
        except error, msg:
            if retries == 0:
                raise error, msg
            if not re.search("locked", str(msg)):
                raise error, msg
            time.sleep(wait)
            retries -= 1
            continue
        break
    return cc

def getColumnNames( dbhandle, table ):
    '''get column names of a table from a database.'''
    
    cc = executewait( dbhandle, "SELECT * FROM %s LIMIT 1" % table )
    return tuple([ x[0] for x in cc.description ])


def getTables( dbhandle ):
    
    cc = executewait( dbhandle, """select name from sqlite_master where type='table'""")
    return tuple([ x[0] for x in cc ])

