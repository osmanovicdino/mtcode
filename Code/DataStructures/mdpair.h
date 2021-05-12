#ifndef MDPAIR_H
#define MDPAIR_H

struct mdpair
{
    int a;
    int b;
    mdpair() : a(0), b(0) {}
    mdpair(const int i) : a(i), b(i) {}
    mdpair(const int &aa, const int &bb) : a(aa), b(bb) {}
    bool operator<(mdpair const &rhs) const
    {
        return a < rhs.a || (a == rhs.a && b < rhs.b);
    }

    mdpair& operator=(const mdpair &rhs) {
        a=rhs.a;
        b=rhs.b;
        return *this;
    }
    
    bool operator==(const mdpair &rhs) const
    {
        if (a == rhs.a && b == rhs.b)
        {
            return true;
        }
        else if (a == rhs.b && b == rhs.a)
        {
            return true;
        }
        else
            return false;
    }

    friend ostream &operator<<(ostream &gt, const mdpair &lhs)
    {
        gt << lhs.a;
        gt << " ";
        gt << lhs.b;
        return gt;
    }
};

bool mdpairCompare(const mdpair &m1, const mdpair &m2) {
    return m1.a < m2.a;
}

struct mdtriplet {
    int a;
    int b;
    int c;

    mdtriplet() : a(0), b(0), c(0) {}
    mdtriplet(const int &aa, const int &bb, const int &cc) : a(aa), b(bb), c(cc) {}
};

struct dispair
{
    int a;
    double b;
    dispair() : a(0), b(0.) {}
    dispair(const int &aa, const double &bb) : a(aa), b(bb) {}
    bool operator<(dispair const &rhs) const
    {
        return b < rhs.b || (b == rhs.b && a < rhs.a);
    }
};

#endif /* MDPAIR_H */
