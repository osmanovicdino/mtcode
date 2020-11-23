#ifndef MDPAIR_H
#define MDPAIR_H

struct mdpair
{
    int a;
    int b;
    mdpair() : a(0), b(0) {}
    mdpair(const int &aa, const int &bb) : a(aa), b(bb) {}
    bool operator<(mdpair const &rhs) const
    {
        return a < rhs.a || (a == rhs.a && b < rhs.b);
    }
    bool operator==(const mdpair &rhs)
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
