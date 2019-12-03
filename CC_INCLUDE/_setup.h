#ifndef __SETUP_H__
#define __SETUP_H__
class _setup {
public: _setup();
public: void init();
public: int getnum_members(char ch);
public: void setxc(int i, double x, char ch);
public: void setyc(int i, double y, char ch);
public: void setzc(int i, double z, char ch);
public: double getxc(int i, char ch);
public: double getyc(int i, char ch);
public: double getzc(int i, char ch);
};
#endif
