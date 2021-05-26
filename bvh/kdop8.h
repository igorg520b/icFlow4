#ifndef KDOP24_H
#define KDOP24_H

namespace icy { class kDOP8; }

class icy::kDOP8
{
public:
    double d[4], g[4]; // 0-3 lower boundaries; 4-7 higher boundaries
    double ctrX, ctrY;
    double dX, dY;

    void Reset();
    bool Overlaps(kDOP8 &b);
    void Expand(double x, double y);
    void Expand(kDOP8 &b);
//    double centerX() { return (d[0] + g[0]) / 2; }
//    double centerY() { return (d[1] + g[1]) / 2; }
//    void Dimensions(double &dx, double &dy)
//    {
//        dx = g[0] - d[0];
//        dy = g[1] - d[1];
//    }

private:
    inline void MinMax(double p, double &mi, double &ma);
};

#endif // KDOP24_H
