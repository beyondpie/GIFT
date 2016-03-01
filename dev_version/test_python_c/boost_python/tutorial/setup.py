from distutils.core import setup, Extension
module1 = Extension('spam', sources = ['test_python_api.c'])
setup (name = 'spam',
       version = '1.0',
       description = 'This is a demo package.',
       ext_modules = [module1] )
