# -*- python -*-

### Python imports
from __future__ import print_function    # python2/3 compatability
import glob

### SCons imports
Import('rootEnv')

def emit_version_header(target, source, env):
    """Populate the version header file with current SVN information."""
    version_h = """
#ifndef _AFFY_VERSION_H_
#define _AFFY_VERSION_H_
#define AFFY_REVISION "%s"
extern const char affy_version[];
#endif
"""
    try:
        import pysvn
        
        c = pysvn.Client()
        subst = c.info(env.Dir('#/').abspath).revision.number
    except ImportError:
        from datetime import datetime
#        subst = '(SVN, rev unknown)'
        now = datetime.now()
        subst = "build date: %s" % now.strftime("%Y-%m-%d %H:%M:%S")
    
    try:
        f = open(str(target[0]), mode='w')
        print (version_h % (subst,), file=f)
        f.close()
    except:
        return 1
    
    return None

### Subdirectories containing libaffy source files.  Each of these will
### be globbed for .c files.
srcsubdirs = ['io', 'mas5', 'rma', 'utils', 'halloc', 'iron_generic', 'iron']
libsrcs = []

for dir in srcsubdirs:
    libsrcs += glob.glob(dir + '/*.c')

### Generate the current version header file.
version_h = rootEnv.Command('include/affy_version.h', 
                            None, 
                            emit_version_header)

#libsrcs.append('swigtest/libaffy.i')

### Build both shared and static versions of libaffy.
rootEnv.StaticLibrary(target = 'affy', source = libsrcs)
#rootEnv.SharedLibrary(target = 'affy', source = libsrcs)

### Documentation target.  Building the docs can be fairly tricky if
### your platform is weird.  So by default leave them alone.
if ARGUMENTS.get('build_docs', 0):
    SConscript('doc/SConscript', exports=['rootEnv'])
