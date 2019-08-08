/* 
 * File:   Strain2D.cpp
 * Author: justin
 * 
 * Created on June 7, 2015, 10:14 PM
 */

#include "Strain2D.h"

namespace ncorr {    
    
// Static factory methods ----------------------------------------------------//
Strain2D Strain2D::load(std::ifstream &is) {
    // Form empty Strain2D then fill in values in accordance to how they are saved
    Strain2D strain;
    
    // Load eyy
    strain.eyy = Data2D::load(is);
	
	// Load exy
    strain.exy = Data2D::load(is);
	
	// Load exx
    strain.exx = Data2D::load(is);
	
	// Load dvy
    strain.dvy = Data2D::load(is);
    
    // Load dvx
    strain.dvx = Data2D::load(is);
    
    // Load duy
    strain.duy = Data2D::load(is);
	
	// Load dux
	strain.dux = Data2D::load(is);
    
    return strain;
}
    
// Operators interface -------------------------------------------------------//
std::ostream& operator<<(std::ostream &os, const Strain2D &strain) {
    os << "Eyy data: " << '\n' << strain.get_eyy();
    os << '\n' << "Exy data: " << '\n' << strain.get_exy();
	os << '\n' << "Exx data: " << '\n' << strain.get_exx();
	os << '\n' << "dvy data: " << '\n' << strain.get_dvx();
	os << '\n' << "dvx data: " << '\n' << strain.get_dvx();
    os << '\n' << "duy data: " << '\n' << strain.get_duy();
	os << '\n' << "dux data: " << '\n' << strain.get_dux();
    
    return os;
}

void imshow(const Strain2D &strain, Strain2D::difference_type delay) { 
    // Just show each separately for now. If you combine into one buffer, you 
    // must scale each as their ranges might be different.
    imshow(strain.eyy, delay);
	imshow(strain.exy, delay);
	imshow(strain.exx, delay);
	imshow(strain.dvy, delay); 
    imshow(strain.dvx, delay); 
    imshow(strain.duy, delay); 
	imshow(strain.dux, delay);
}  

bool isequal(const Strain2D &strain1, const Strain2D &strain2) {
    return isequal(strain1.eyy, strain2.eyy) && isequal(strain1.exy, strain2.exy) && isequal(strain1.exx, strain2.exx) && isequal(strain1.dvy, strain2.dvy) && isequal(strain1.dvx, strain2.dvx) && isequal(strain1.duy, strain2.duy) && isequal(strain1.dux, strain2.dux);
}

void save(const Strain2D &strain, std::ofstream &os) {        
    // Save eyy -> exy -> exx
    save(strain.eyy, os);
	save(strain.exy, os);
	save(strain.exx, os);
	save(strain.dvy, os);
    save(strain.dvx, os);
    save(strain.duy, os);
	save(strain.dux, os);
}

// Interpolator --------------------------------------------------------------//
Strain2D::nlinfo_interpolator Strain2D::get_nlinfo_interpolator(difference_type region_idx, INTERP interp_type) const {
    return details::Strain2D_nlinfo_interpolator(*this, region_idx, interp_type);
}
    
}
