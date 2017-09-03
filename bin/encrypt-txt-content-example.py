#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys
import getpass

validPassword = 'password'

inputPassword = getpass.getpass(prompt="Key: ")

if inputPassword == validPassword:
    print 'You have access!'
    print 'Welcome!'
    mybool = raw_input('Write contents into txt?[y/n]')
    if mybool == 'y':
    	secret_content = u""" SECRET CONTENT """
        f = open('hstory01.txt', 'w')
        f.write(secret_content.encode('UTF-8'))
        f.close()
    else:
        print "残念~ 下次吧！"
else:
    print 'Access denied!'
    sys.exit(0)

