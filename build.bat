@echo off

WHERE cl 
IF %ERRORLEVEL% NEQ 0 call %VCVARSALL% x64

set DEBUG=1

setlocal ENABLEDELAYEDEXPANSION

set WIGNORE=-wd4201 -wd4505 -wd4996 -wd4100 -wd4238

set CommonCompilerFlags=/MD -fp:fast -fp:except- -nologo /O2 -Oi -W4 -Gm- -GR- -EHa -FC -Z7 /FeRTTut -DDEBUG=%DEBUG% %WIGNORE% 
set CommonLinkerFlags=shell32.lib 
set ExtraLinkerFlags=/NODEFAULTLIB:"LIBCMT" -incremental:no -opt:ref /ignore:4099

IF %DEBUG% NEQ 1 set Console=/SUBSYSTEM:windows /ENTRY:mainCRTStartup

IF NOT EXIST bin mkdir bin
pushd bin

del *.pdb > NUL 2> NUL

echo Compilation started  on: %time%
cl %CommonCompilerFlags% ..\src\main.cpp /link %ExtraLinkerFlags% %CommonLinkerFlags% %Console%

echo Compilation finished on: %time%
popd
