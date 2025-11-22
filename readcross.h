/**
 * readcross.h
 * 
 * Header file for gas opacity reading and reinterpolation functions
 */

#ifndef READCROSS_H
#define READCROSS_H

// Opacity reinterpolation functions
extern void reinterpolate_all_opacities();
extern void cleanup_opacity_cache(void);

// Main opacity reading function
void readcross(char Fname[], double **xsc);
void read_all_opacities();

#endif // READCROSS_H

