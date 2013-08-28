
Windows* OS Installation Notes:
-------------------------------
    
    Extract all files in iaca-win32.zip into a directory on your machine
    
    Add the complete path <full path of your directory >\iaca-win32 to the PATH environment variable
    For example:
    set PATH = %PATH%;"C:\iaca-win32"

Linux* OS Installation Notes:
-----------------------------

    Extract all files in iaca-lin32.zip into a directory on your machine

    Add the complete path <full path of your directory >/iaca-lin32/bin to the PATH environment variable
    For example:
    set PATH = $PATH:/home/iaca-lin32/bin"

    Add the complete path <full path of your directory >/iaca-lin32/lib to the LD_LIBRARY_PATH environment variable
    For example:
    set LD_LIBRARY_PATH = $LD_LIBRARY_PATH:/home/iaca-lin32/lib"

Mac* OS Installation Notes:
---------------------------

    Extract all files in iaca-mac32.zip into a directory on your machine

    Add the complete path <full path of your directory >/iaca-mac32/bin to the PATH environment variable
    For example:
    set PATH = $PATH:/home/iaca-mac32/bin"

    Add the complete path <full path of your directory >/iaca-mac32/lib to the DYLD_LIBRARY_PATH environment variable
    For example:
    set DYLD_LIBRARY_PATH = $DYLD_LIBRARY_PATH:/home/iaca-mac32/lib"

For further details on operating Intel(R) Architecture Code Analyzer please refer 
to the user's guide that is part of the release.

Intel(R) Architecture Code Analyzer 2.0.1 Release Notes
-----------------------------------------------------------

- Fixed a bug where –graph option failed to produce graph file.
