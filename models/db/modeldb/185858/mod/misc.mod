: $Id: misc.mod,v 1.24 2011/10/14 15:00:26 samn Exp $

COMMENT
Misc. routines:
sassign() // assign a string from system
dassign()// assign a double
nokill() // chatch SIGHUP
prtime() // gives date/time
fspitchar(c,file) // sends single char to a file
spitchar(c)       // sends single char to stdout: eg c=1 => ^A
file_exist(file) // returns 1 if filename exists
hocgetc(file) // get single char from a file

  Note that with a SUFFIX equal to "nothing" these functions do not
have a suffix in hoc.  Thus to call sassign() in hoc use simply type
"sassign()" <- without the quotes.

    file_exist(filename)
        - returns 1 if filename exists

    sassign()  (string assign, written by Bill Lytton)
        - This routine is used to set a string in Hoc to something that has
          been returned by a system call.  sassign("name","shell_call ...")
          will produce a file called "sassign" in the cwd that will contain
          a hoc call that sets string 'name' to the result of shell_call 
          which should be a string.
        
    dassign()  (double assign, written and used by Bill Lytton)
        - This routine is used to set a variable in Hoc to something that has
          been returned by a system call.  sassign("name","shell_call ...")
          will produce a file called "dassign" in the cwd that will contain
          a hoc call that sets variable 'name' to the result of shell_call 
          which should be a number.

ENDCOMMENT
                           
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX nothing
}

FUNCTION istmpobj () {
VERBATIM
  _listmpobj=hoc_is_tempobj_arg(1);
ENDVERBATIM  
}

