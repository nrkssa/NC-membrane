#ifndef __DEBUGGER_H__
#define __DEBUGGER_H__
class _debugger{
public:
	stringstream emailid,st1,bash_cmd;

public: _debugger(){}
public: ~_debugger(){}
public: void init();
public: void write_mailscript(string subject);
public: void send_email(string subject);
};
#endif
