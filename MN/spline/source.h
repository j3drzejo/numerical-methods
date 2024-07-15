#ifndef SPLINE_H
#define SPLINE_H

class spline {
public:
    spline(int n);
    ~spline();
    void set_points(double x[], double y[]);
    double operator()(double z) const;

private:
    int n;
    double *x;
    double *y;
    double *h;
    double *a;
    double *b;
    double *c;
    double *d;
};

#endif // SPLINE_H
