# -*- python -*-

### Python imports
import platform

### Constant data structures

## Adjust special compiler flags here.  Each compiler tool has its own
## bundles of named flags.  However, the idea is that each should
## support a few common bundles (optimize, debug, etc.) alongside
## compiler-specific bundles.  Bundle selection is done on the command
## line.
cflags_dict = { 'intelc': { 'optimize': ['-O3', '-ipo', '-xN'],
                            'debug': ['-g'], 
                            'strict': ['-Wall', '-strict-ansi'], },
                'gcc': { 'strict_gcc2x': ['--static', '-ansi', '-pedantic', '-Wall'],
                         'strict_gcc3x': ['--static', '-std=c99', '-pedantic', '-Wall'],
                         'strict': ['--static', '-pedantic', '-Wall'],
                         'optimize': ['--static', '-O3', '-finline', '-funroll-loops', '-pedantic', '-Wall', '-Wno-unused'],
                         'datextractor': ['-DSTORE_CEL_QC', '--static', '-O3', '-finline', '-funroll-loops'],
                         'profile': ['--static', '-g', '-pg'],
                         'debug': ['-g'],
                         'debugo3': ['-g', '-O3', '-fno-omit-frame-pointer'], },
                'cygmingw': { 'strict_gcc2x': ['--static', '-ansi', '-pedantic', '-Wall'],
                              'strict_gcc3x': ['--static', '-std=c99', '-pedantic', '-Wall'],
                              'strict': ['--static', '-pedantic', '-Wall'],
                              'optimize': ['--static', '-O3', '-finline', '-funroll-loops', '-pedantic', '-Wall', '-Wno-unused'],
                              'datextractor': ['-DSTORE_CEL_QC', '--static', '-O3', '-finline', '-funroll-loops'],
                              'debug': ['--static', '-g'], },
                'msvc': { 'optimize': ['/Ox'],
                          'debug': ['/Od', '/Zi'],
                          'strict': ['/Wall', '/Zi'], },
                }

### Help text
help_text = """
Supported values for `compiler' and `cflags_bundle':\n
"""

#for tool, d in cflags_dict.iteritems():
#
# items() is inefficient, but at least it is compatible with both 2/3
#  we don't care about efficiency here anyways
for (tool, d) in cflags_dict.items():
    help_text = help_text + "%15s: %s\n" % (tool,
                                            ' '.join(d.keys()))

Help(help_text)

### Check for specified compiler tool or cflags bundle on the command line.
alt_compiler = ARGUMENTS.get('compiler', None)
cflags_bundle = ARGUMENTS.get('cflags_bundle', None)

tool_list = ['default']

if alt_compiler:
    tool_list.append(alt_compiler)

rootEnv = Environment(tools=tool_list)

### Preprocessor defines which are compiler-independent
cpp_defines = \
    {
    platform.system(): None,
## Uncomment below to remove assertions
#    'NDEBUG': None
    }


### Compiler configuration section

## Simple-minded way of guessing the tool package based on the value
## of 'CC' in the environment.  Will NOT always work properly.  This
## is currently kind of a weakness in SCons.
def guess_compiler_tool(cc_str, platform):
    if cc_str == 'gcc':
        if platform == 'cygwin':
            return 'cygmingw'
        else:
            return 'gcc'
    elif cc_str == 'xlc':
        return 'aixcc'
    elif cc_str == 'bcc32':
        return 'bcc32'
    elif cc_str == 'icc' or cc_str == 'icl':
        return 'intelc'
    elif cc_str == 'mwcc':
        return 'mwcc'
    elif cc_str == 'cl':
        return 'msvc'


if cflags_bundle:
    if alt_compiler:
        compiler_name = alt_compiler
    else:
        compiler_name = guess_compiler_tool(rootEnv['CC'],
                                            rootEnv['PLATFORM'])
    rootEnv.Append(CCFLAGS = cflags_dict[compiler_name][cflags_bundle])


# Compile statically linked gcc binaries, for maximum binary portability.
# To avoid various problems, set CC rather than adding LINKFLAGS.
# Add -w to cygwin to get rid of "warning: ISO C requires whitespace after the macro name"
#  I couldn't figure out which -Wno- option to toggle, so I just kill them all
if cflags_bundle != 'debug':
  if cflags_bundle != 'profile':
    if rootEnv['CC'] == 'gcc':
      rootEnv['CC'] = 'gcc'
    if rootEnv['PLATFORM'] == 'cygwin':
      rootEnv['CC'] = 'i686-pc-mingw32-gcc.exe -w'
  if cflags_bundle == 'profile':
    if rootEnv['CC'] == 'gcc':
      rootEnv['CC'] = 'gcc -pg'
    if rootEnv['PLATFORM'] == 'cygwin':
      rootEnv['CC'] = 'i686-pc-mingw32-gcc.exe -pg -w'


### Configuration checks

# Rudimentary check for POSIX system
confCtx = Configure(rootEnv)

# check for POSIX first, since cygwin is more POSIX than WIN32
if confCtx.CheckCHeader('unistd.h'):
    print ("You appear to be compiling on a POSIX system.")
    cpp_defines['AFFY_POSIX_ENV'] = None
elif (confCtx.CheckCHeader('windows.h') and
    confCtx.CheckDeclaration('WIN32', includes='#include <windows.h>')):
    print ("You appear to be compiling on a Win32 system.")
    cpp_defines['AFFY_WIN32_ENV'] = None
else:
    print ("Couldn't detect a supported OS environment, build will likely fail")
if confCtx.CheckCHeader('netcdf.h'):
    print ("NetCDF package seems to be available.")
    cpp_defines['AFFY_HAVE_NETCDF'] = None

print ("DEBUG %s %s" % (confCtx.CheckCHeader('unistd.h'), rootEnv['PLATFORM']))

confCtx.Finish()


rootEnv.Replace(CPPPATH = ['#/libutils/', 
                           '#/libutils/argp/', 
                           '#/libutils/getopt/', 
                           '#/libaffy/include/',
                           '#/libaffy/halloc/'],
                LIBPATH = ['#/libutils/', '#/libaffy/'],
                CPPDEFINES = cpp_defines)


SConscript(dirs = ['libaffy', 'affy-apps'], exports=['rootEnv'])


### Import the libutils compilation environment and modify it to be
### built subsidiary to this SConstruct.

utilsEnv = SConscript('libutils/SConstruct')

if cflags_bundle:
    utilsEnv.Replace(CCFLAGS = cflags_dict[compiler_name][cflags_bundle])

utilsEnv.Replace(CPPPATH = ['#/libutils/', 
                            '#/libutils/argp', 
                            '#/libutils/getopt'],
                 CPPDEFINES = cpp_defines)
