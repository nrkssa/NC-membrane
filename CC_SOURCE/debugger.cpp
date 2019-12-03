#include "_declarations.h"
#include "_debugger.h"

void _debugger :: init(){
	emailid<<"ramn@seas.upenn.edu";
	st1<<"./mailerscript.sh";
	bash_cmd<<"bash "<<st1.str().c_str();
}

void _debugger  :: write_mailscript(string subject){
	ofstream so1;
	so1.open(st1.str().c_str(),ios::out);
	so1<<"#!/bin/bash"<<endl;
	so1<<"str=`pwd`"<<endl;
	so1<<"echo \"Error message: $str\""<<" "<<subject <<
	     "| mail -s \"Error message from location ($str)\" "<<emailid.str().c_str()<<endl;         //writes a bash script with contents of the mail.
	so1.close();
}

void _debugger  :: send_email(string subject){                               // send an email to emailid with the subject passed from the code.
	write_mailscript(subject);
	system(bash_cmd.str().c_str());
}

