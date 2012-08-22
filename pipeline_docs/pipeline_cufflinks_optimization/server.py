#!/bin/env python

import web, os

from SphinxReport import Cache
from SphinxReport import Utils
from SphinxReport import DataTree

urls = ( '/data/(.*)', 'DataTable',
         '/index/(.*)', 'Index'  )

# expose zip within templates
global_vars = {'zip': zip }

render = web.template.render('templates/', globals = global_vars)

app = web.application(urls, globals() )

class DataTable:
    '''render data retrieved from cache as a table.'''

    def GET(self, tracker):

        cache = Cache.Cache( tracker, mode = "r" )
        data = DataTree.fromCache( cache )
        table, row_headers, col_headers = DataTree.tree2table( data )

        return render.data_table(table, row_headers, col_headers )


if __name__ == "__main__": app.run()



