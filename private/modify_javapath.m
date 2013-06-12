function modify_javapath()
    javaaddpath(strrep(mfilename('fullpath'), mfilename(), ''))
