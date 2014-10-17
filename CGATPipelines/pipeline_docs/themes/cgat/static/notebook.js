/* 
   notebook.js

   Author: Andreas Heger
   7/10/2014

   Instructions:
   Download this file
   Add <script src="notebook.js"></script> to your HTML
   Add class="sortable" to any table you'd like to make sortable
   Click on the headers to sort

   Licenced as X11: http://www.kryogenix.org/code/browser/licence.html
   This basically means: do what you want with it.
*/

function setCookie(cname, cvalue, exdays) {
    var d = new Date();
    d.setTime(d.getTime() + (exdays*24*60*60*1000));
    var expires = "expires="+d.toUTCString();
    document.cookie = cname + "=" + cvalue + "; " + expires;
}

function getCookie(cname) {
    var name = cname + "=";
    var ca = document.cookie.split(';');
    for(var i=0; i<ca.length; i++) {
        var c = ca[i];
        while (c.charAt(0)==' ') c = c.substring(1);
        if (c.indexOf(name) != -1) return c.substring(name.length, c.length);
    }
    return "";
}

function checkCookie() {
    var url = getCookie("NotebookURL");
    if (url == "") {
	url = "http://localhost:8888";
    }
    url = prompt("Please enter URL:", url);
    if (url != "" && url != null) {
        setCookie("NotebookURL", url, 365);
    }
    return url
}


function create_notebook(option_string)
{
    url = checkCookie()

    if (url == null) { return };

    var data = {
	content: {
	    metadata: {
		name:""
            },
	    nbformat: 3,
	    nbformat_minor: 0,
	    "worksheets": [
		{
		    "cells": [
			{
			    "cell_type": "code",
			    "collapsed": false,
			    "input": [
				'%matplotlib inline',
				'import CGATReport.test',
				'result = CGATReport.test.main('+option_string+')',
			    ],
			    "language": "python",
			    "metadata": {},
			    "outputs": [
				{
				    "output_type": "stream",
				    "stream": "stdout",
				    "text": []
				}
			    ],
			    "prompt_number": 1
			},
		    ],
		    "metadata": {}
		}
	    ]
	}
    };

    var msg = JSON.stringify(data);
    url = url + "/api/notebooks"
    // alert("I am POSTing this to " + url + ":\n\n" + msg);
    $.post(url, msg);
}

