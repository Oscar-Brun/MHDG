#!/bin/bash
# --- Chemins vers tes installations locales ---
export MHDG_GMSH_DIR=/home/obrun/logiciels/gmsh-4.11.1-Linux64-sdk
export MHDG_PASTIX_DIR=/home/obrun/logiciels/pastix_install
export MHDG_SCOTCH_DIR=/usr

# --- Configuration pour PaStiX ---
export PKG_CONFIG_PATH=$MHDG_PASTIX_DIR/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$MHDG_PASTIX_DIR/lib:$MHDG_GMSH_DIR/lib:$LD_LIBRARY_PATH

echo "Environnement MHDG chargé pour WEST !"


