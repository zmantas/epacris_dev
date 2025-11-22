/**
 * readcia.h
 * 
 * Header file for Collision-Induced Absorption (CIA) opacity functions
 */

#ifndef READCIA_H
#define READCIA_H

// CIA opacity reinterpolation and cleanup functions
extern void reinterpolate_all_cia_opacities();
extern void cleanup_cia_cache();

#endif // READCIA_H

