#include "molstate.h"

/*! \file molstate.C
\brief Molecular state: Stores Geometry, Normal Modes & Frequencies. 
Also reads input data from the XML file (vm'06) 
\ingroup DATA_CLASSES
*/

//------------------------------
MolState::MolState () {
  
  momentOfInertiaTensor.Adjust(CARTDIM,CARTDIM);
  if_aligned_manually=false;
  if_nm_reordered_manually=false;
  normModesOrder.clear();
}


//------------------------------
MolState::MolState (const MolState& other) {

  atoms=other.atoms;
  normModes=other.normModes;
  gradient=other.gradient;
  ifLinear=other.ifLinear;
  energy=other.energy;
  centerOfMass = other.centerOfMass;
  momentOfInertiaTensor = other.momentOfInertiaTensor;
  if_aligned_manually = other.if_aligned_manually;
  if_nm_reordered_manually =other.if_nm_reordered_manually;
  normModesOrder = other.normModesOrder;
  reduced_masses=other.reduced_masses;
}

//------------------------------
MolState& MolState::operator=(const MolState& other) {
  if (this != &other) {
    
    atoms=other.atoms;
    normModes=other.normModes;
    gradient=other.gradient;
    ifLinear=other.ifLinear;
    energy=other.energy;
    centerOfMass = other.centerOfMass;
    momentOfInertiaTensor = other.momentOfInertiaTensor;
    if_aligned_manually = other.if_aligned_manually;
    if_nm_reordered_manually =other.if_nm_reordered_manually;
    normModesOrder = other.normModesOrder;
    reduced_masses=other.reduced_masses;
  }
  return *this;
}


//------------------------------
void MolState::align()
{
  // Shift center of mass to the origin:
  shiftCoordinates( getCenterOfMass() );
  /* const char * center_of_mass_name = "Center of mass: "; */
  /* getCenterOfMass().print(center_of_mass_name); */


  // Aligh moment of inertia principal axes:
  KMatrix MOI_tensor(3,3), MOI_eigenValues(3,1), MOI_eigenVectors(3,3);
  MOI_tensor = getMomentOfInertiaTensor();
  //MOI_tensor.Print("Moment of inertia tensor");
  MOI_eigenVectors=MOI_tensor.Herm_Diag(MOI_eigenValues, true); // true=>eigen vals in the descending order; eigen vectors in ROWS.
  // MOI_eigenVectors.Print("MOI eigen vectors");
  // MOI_eigenValues.Print("MOI eigen values");
  
  // compute the determinant of MOI matrix:
  double det_tmp =  MOI_eigenVectors.Determinant();
  // if determinant is = 1 than it is a proper rotation; 
  // if it is = -1 than it is rotation+reflection (switch from left handed to the right handed coord.system); swap x&y axes:
  if (det_tmp < 0)
    MOI_eigenVectors.SwapRows(0,1);

  transformCoordinates(MOI_eigenVectors); // rotate coordinates using matrix of the MOI eigen vectors 
  std::cout << "Coordinate axes were aligned with MOI principal axes; center mass was shifted to the origin.\n";
}


//------------------------------
void MolState::align(MolState& other)
{
  //align the state with the "other" state by rotating around every axes (x,y,z) by pi/2;
  // by minimizing the "sum" of distances between the same atoms: SUM[i=1..Natoms](deltaRi)
  // DeltaRi: distance between i'th atoms of the initial and target states

  double min_diff=DBL_MAX; // minimum of the "sum"
  int min_x_rot, min_y_rot, min_z_rot; // rotations, that correspond to the minimum of the "sum"

  for (int z_rot=0;  z_rot<4; z_rot++)
  {
    for (int y_rot=0;  y_rot<4; y_rot++)
    {
      for (int x_rot=0;  x_rot<4; x_rot++)
      {
        double diff=getGeomDifference(other);
        if (diff<=min_diff)
        {
          min_diff=diff;
          min_x_rot=x_rot;
          min_y_rot=y_rot;
          min_z_rot=z_rot;
        }
        rotateX_90deg(); // after 4 rotations by pi/2, it gets the "0" rotation geometry and starts over.
      }
      rotateY_90deg(); // after 4 rotations by pi/2, it gets the "0" rotation geometry and starts over.
    }
    rotateZ_90deg(); // after 4 rotations by pi/2, it gets the "0" rotation geometry.
  }// at the end the molecule rotated by 4*pi/2 around x and y, i.e. it is non-rotated
  
  // rotate to the angles with the min norm of the deltaRi differences (min_x_rot, min_y_rot, min_z_rot):
  rotate(min_x_rot*PI/2.0, min_y_rot*PI/2.0, min_z_rot*PI/2.0);

  std::cout << "Also rotated by " << min_z_rot <<"/2*pi, " 
    << min_y_rot <<"/2*pi, and " << min_x_rot <<"/2*pi CCW around z, y, and x.\n";
  std::cout << "The norm of the geometry difference from the initial state is " << sqrt(min_diff)<<"\n";
}


//------------------------------
bool MolState::ifSimilar(MolState& other)
{
  bool return_bool=true;
  std::string error_str="";
  if ( NAtoms()!=other.NAtoms() )
    {
      std::cout << "\nError: different number of atoms in the initial and the target states\n\n";
      return_bool=false;	
    }
  else if ( NNormModes()!=other.NNormModes() )
    {
      std::cout << "\nError: different number of normal modes in the initial and the target state\n\n";
      return_bool=false;	
    }

  return return_bool;
}


//------------------------------
Vector3D& MolState::getCenterOfMass()
{
  centerOfMass.reset();
  double totalMass=0;
  for (int i=0; i<atoms.size();i++)
    {
      for (int axis=0; axis<3; axis++)
	centerOfMass.getCoord(axis)+=atoms[i].getMomentumProj(axis);
      totalMass+=atoms[i].getMass();
    }
  centerOfMass*=1/totalMass;
  return centerOfMass;
}

//------------------------------
KMatrix& MolState::getMomentOfInertiaTensor()
{
  //  I_ij=SUM_k[m_k*(r_k^2*delta_ij-r_ki*r_kj)]
  for (int i=0; i< CARTDIM; i++)
    for (int j=0; j<CARTDIM; j++)
      {
	momentOfInertiaTensor.Elem2(i,j)=0;
	for (int atom=0; atom<NAtoms(); atom++)
	  { 
	    momentOfInertiaTensor.Elem2(i,j)-=atoms[atom].getCoord(i)*atoms[atom].getCoord(j)*atoms[atom].getMass();

	    // add R^2, i.e. I_zz=m*(R^2-r_z^2)=m*(r_x^2+r_y^2)
	    if (i==j)
	      momentOfInertiaTensor.Elem2(i,j)+=atoms[atom].getR()*atoms[atom].getR()*atoms[atom].getMass();
	  }
      }
  return momentOfInertiaTensor;
}


//------------------------------
void MolState::shiftCoordinates(Vector3D& vector)
{
  for (int i=0; i<atoms.size(); i++)
    atoms[i].shiftCoordinates(vector);
  // there is no need to shift normal coordinates!!
}

//------------------------------
void MolState::transformCoordinates(const KMatrix& matrix_3x3)
{
  for (int i=0; i<atoms.size(); i++)
    atoms[i].transformCoordinates(matrix_3x3);
  for (int i=0; i<normModes.size(); i++)
    normModes[i].transformCoordinates(matrix_3x3);
}

//------------------------------
void MolState::rotateX_90deg()
{
  for (int i=0; i<atoms.size(); i++)
    atoms[i].rotateX_90deg();
  for (int i=0; i<normModes.size(); i++)
    normModes[i].rotateX_90deg();
}

//------------------------------
void MolState::rotateY_90deg()
{
  for (int i=0; i<atoms.size(); i++)
    atoms[i].rotateY_90deg();
  for (int i=0; i<normModes.size(); i++)
    normModes[i].rotateY_90deg();
}

//------------------------------
void MolState::rotateZ_90deg()
{
  for (int i=0; i<atoms.size(); i++)
    atoms[i].rotateZ_90deg();
  for (int i=0; i<normModes.size(); i++)
    normModes[i].rotateZ_90deg();
}

//------------------------------
void MolState::rotate(const double alpha_x, const double alpha_y, const double alpha_z)
{
  // three rotation matrices (instead of making one matrix) around x, y, z axes
  KMatrix R(CARTDIM,CARTDIM), Rx(CARTDIM,CARTDIM), Ry(CARTDIM,CARTDIM), Rz(CARTDIM,CARTDIM);

  Rx.Set(0);
  Rx.Elem2(0,0)=1;
  Rx.Elem2(1,1)=cos(alpha_x);
  Rx.Elem2(2,2)=cos(alpha_x);
  Rx.Elem2(2,1)=-sin(alpha_x);
  Rx.Elem2(1,2)=sin(alpha_x);

  Ry.Set(0);
  Ry.Elem2(1,1)=1;
  Ry.Elem2(0,0)=cos(alpha_y);
  Ry.Elem2(2,2)=cos(alpha_y);
  Ry.Elem2(2,0)=sin(alpha_y);
  Ry.Elem2(0,2)=-sin(alpha_y);

  Rz.Set(0);
  Rz.Elem2(2,2)=1;
  Rz.Elem2(0,0)=cos(alpha_z);
  Rz.Elem2(1,1)=cos(alpha_z);
  Rz.Elem2(1,0)=-sin(alpha_z);
  Rz.Elem2(0,1)=sin(alpha_z);

  // Rotations are applied in the order: x, y, and z
  // which gives the overall rotation matrix
  // R = Rz Ry Rx
  R.SetDiagonal(1);
  R*=Rx;
  R*=Ry;
  R*=Rz;
  
  // and now rotates using matix R:
  transformCoordinates(R);
}
//------------------------------
bool MolState::ifAlignedManually()
{
  if (if_aligned_manually)
    return true;
  else
    return false;
}
//------------------------------
bool MolState::ifNMReorderedManually()
{
  if (if_nm_reordered_manually)
    return true;
  else
    return false;
}

//------------------------------
double MolState::getGeomDifference(MolState& other)
{
  double diff=0;
  for (int i=0; i<atoms.size(); i++)
    for (int j=0; j<CARTDIM; j++) 
      // diff+=DeltaR^2[=DetaX^2+DeltaY^2+DeltaZ^2]
      diff+= (getAtom(i).getCoord(j)-other.getAtom(i).getCoord(j)) * (getAtom(i).getCoord(j)-other.getAtom(i).getCoord(j));
  return diff;
}

//------------------------------
void MolState::applyCoordinateThreshold(const double threshold)
{
  for (int i=0; i<atoms.size(); i++)
    atoms[i].applyCoordinateThreshold(threshold);
  for (int i=0; i<normModes.size(); i++)
    normModes[i].applyCoordinateThreshold(threshold);
}

bool MolState::getNormalModeOverlapWithOtherState(MolState& other, KMatrix& overlap, std::vector<int>& normal_modes_list)
{
// SG: next line not used.
//  overlap; //#NM x #NM
  overlap.Adjust(NNormModes(),NNormModes());
  overlap.Set(0.0);

  //normalization of initial and target normal modes (stored as un-mass-weighted)
  double norm_ini, norm_targ;
  // norm of the L
  for (int nm1=0; nm1<NNormModes(); nm1++)
    for (int nm2=0; nm2<NNormModes(); nm2++)
      {
	norm_ini = 0;
	norm_targ = 0;
	for (int a=0; a<NAtoms(); a++)
	  for (int i=0; i<CARTDIM; i++ )
	    {
	      //x1*x1+y1*y1+..
	      norm_ini+= getNormMode(nm1).getDisplacement().Elem1(a*CARTDIM+i) * getNormMode(nm1).getDisplacement().Elem1(a*CARTDIM+i);
	      //x2*x2+y2*y2+..
	      norm_targ+= other.getNormMode(nm2).getDisplacement().Elem1(a*CARTDIM+i) * other.getNormMode(nm2).getDisplacement().Elem1(a*CARTDIM+i);
	      //x1*x2+y1*y2+...
	      overlap.Elem2(nm1,nm2)+= getNormMode(nm1).getDisplacement().Elem1(a*CARTDIM+i) * other.getNormMode(nm2).getDisplacement().Elem1(a*CARTDIM+i);
	    }
	overlap.Elem2(nm1,nm2)/= sqrt(norm_ini) * sqrt(norm_targ);
      }
  bool return_bool=true;


  // scan the diagonal; if the diagonal element is NOT the maximum element in the raw and the colum -- than ADD this nm to the list;
  double max;
  for (int nm=0; nm<NNormModes(); nm++)
    {
      max=fabs(overlap.Elem2(nm,nm)); // maximum should be at the diagonal
      for (int nm_scan=0; nm_scan<NNormModes(); nm_scan++) // scan row&column
	if ( (fabs(overlap.Elem2(nm,nm_scan))>max) or (fabs(overlap.Elem2(nm_scan,nm))>max) ) // check nm-th row & column
	  {
	    return_bool=false;
	    normal_modes_list.push_back(nm);
	    normal_modes_list.push_back(nm_scan); // the column/row with larger element should be also included
	  }
    }

  // Now Sort and Remove duplicates from the normal_modes_list:
  std::sort( normal_modes_list.begin(), normal_modes_list.end() );
  std::vector<int>::iterator new_end_pos;
  new_end_pos = std::unique( normal_modes_list.begin(), normal_modes_list.end() );
  normal_modes_list.erase( new_end_pos, normal_modes_list.end() );

  return return_bool;
 }



void MolState::Print()
{
  std::cout << "Geometry=\n";
  printGeometry();    
  std::cout << "Normal modes=\n";
  for (int k=0; k<NNormModes(); k++)
    {
      std::cout << "   Frequency="; 
      std::cout << getNormMode(k).getFreq() << '\n';
      std::cout << "   Displacement=\n";
      for (int i=0; i<NAtoms(); i++)
	{
	  for (int l=0; l<CARTDIM; l++)
	    std::cout << getNormMode(k).getDisplacement()[i*CARTDIM+l] << ' ';
	  std::cout << '\n';
	}
    }
  std::cout <<" end of the electronic state \n";   
}

void MolState::printGeometry()
{
  for (int i=0; i<NAtoms(); i++)
    {
      std::cout << std::setw(4) << std::right  << getAtom(i).Name();
      for (int k=0; k<CARTDIM; k++)
	std::cout << std::setw(12) << std::right << std::fixed << std::setprecision(4) << std::showpoint << getAtom(i).Coord(k) << ' '; 
      std::cout << '\n';
    }
}


void MolState::printNormalModes()
{
  // print in qchem format (3 per line)
  int nModesPerLine=3;  // three number of vib. modes per Line
  int nLines = NNormModes() / nModesPerLine;
  if ( NNormModes() % nModesPerLine != 0)
    nLines++;

  for (int n = 0; n < nLines; n++)  // number of blocks with 3 norm.modes. ("lines")
    {
      int current_nModesPerLine = nModesPerLine;
      // for the last entree, nModesPerString may differ from 3
      if (nLines - 1 == n)
  	if ( NNormModes() % nModesPerLine != 0 )
	  current_nModesPerLine = NNormModes() % nModesPerLine;
 
      for (int a=0; a < NAtoms(); a++)   
	{
	  for (int j=0; j < current_nModesPerLine; j++) 
	    {
	      for (int k=0; k < CARTDIM; k++)
	      
		std::cout << std::setw(7) << std::right << std::fixed << std::setprecision(3)
			  << getNormMode(n*nModesPerLine+j).getDisplacement()[a*CARTDIM+k]*sqrt( reduced_masses[n*nModesPerLine+j]/getAtom(a).Mass() )<<' ';
	      std::cout <<  "  ";
	    }
	  std::cout << "\n";
	}
      std::cout << "\n";
    }

  /*

  for (int nm=0; nm < NNormModes(); nm++)
    {
      for (int a=0; a < NAtoms(); a++) 
	{
	  for (int k=0; k < CARTDIM; k++)
	    std::cout << std::setw(10) << std::right << std::fixed << std::setprecision(4)
		      << getNormMode(nm).getDisplacement()[a*CARTDIM+k]*sqrt( reduced_masses[nm]/getAtom(a).Mass() )<<' ';
	  std::cout << '\n';
	}
      std::cout << '\n';
    }

  */


}

void MolState::printGradient()
{
  for (int i=0; i<NAtoms(); i++)
    {
      std::cout << std::setw(4) << std::right  << getAtom(i).Name();
      for (int k=0; k<CARTDIM; k++)
	std::cout << std::setw(12) << std::right << std::fixed << std::setprecision(4) << std::showpoint << gradient.Elem2(i * CARTDIM + k, 0) << ' '; 
      std::cout << '\n';
    }
}



//------------------------------
bool MolState::ifLetterOrNumber(char Ch)
{
if ( ((int(Ch)<=int('Z'))&&(int(Ch)>=int('A')))  
  || ((int(Ch)<=int('z'))&&(int(Ch)>=int('a'))) 
  || ((int(Ch)<=int('9'))&&(int(Ch)>=int('0'))) )
  return true;
else return false;
}


//------------------------------ 
/* Only this function needs to deal with input: 
 -  a node pointing out to initial_state or target_state section in the input
 -  a file where all masses are tabulated
*/
bool MolState::Read(xml_node& node_state, xml_node& node_amu_table)
{
  int i,j,k,l;

  //------------ Read IP (if provided) ----------------------------------
  energy = 0.0; //units are eV
  if(node_state.find_subnode("excitation_energy")) {
    std::string energy_text = "Adiabatic excitation energy = ";
    bool gradient_is_available = node_state.find_subnode("gradient");
    if (gradient_is_available) {
      energy_text = "Vertical excitation energy = ";
    }

    xml_node node_ip(node_state, "excitation_energy", 0);
    std::string units=node_ip.read_string_value("units");
    energy=node_ip.read_node_double_value();
    std::cout << std::fixed << std::setprecision(6); //  <<  std::setw(6);
    std::cout << energy_text << energy << " " << units << std::endl;
    
    if ( !covert_energy_to_eV(energy,units) ) {
      std::cout << "\nError! Unknow units of the IP: \"" << units <<"\"\n  (should be equal to \"eV\", \"K\", or \"cm-1\")\n\n";
      exit(1);
    }

    std::cout << energy_text << energy << " eV " << std::endl;
  }
    

  //------------ Read the Geometry --------------------------------------
  xml_node node_geometry(node_state,"geometry",0);
  
  int tmp_nAtoms, tmp_nNormMds;
  tmp_nAtoms=node_geometry.read_int_value("number_of_atoms");
  
  std::string units=node_geometry.read_string_value("units");
  
  ifLinear= node_geometry.read_bool_value("linear"); 
  if (ifLinear) 
    tmp_nNormMds = (3*tmp_nAtoms - 5);
  else
    tmp_nNormMds = (3*tmp_nAtoms - 6);

  My_istringstream geom_iStr(node_geometry.read_string_value("text"));
  atoms.clear();
  
  Atom tmp_atom;
  double coeff=(units=="au") ? AU2ANGSTROM : 1.0;
  for (i=0; i<tmp_nAtoms; i++) {

    std::string tmp_atomName;
    //Get atomic name:
    geom_iStr.getNextWord(tmp_atomName);
    tmp_atom.Name() = tmp_atomName;

    //get coordinates & convert to a.u.
    // TODO: `coeff` converts to Angstroms. The above comment is confusing. Pawel Feb '22
    tmp_atom.Coord(0)=geom_iStr.getNextDouble()*coeff;
    tmp_atom.Coord(1)=geom_iStr.getNextDouble()*coeff;
    tmp_atom.Coord(2)=geom_iStr.getNextDouble()*coeff;
    
    if(geom_iStr.fail()) {
      std::cout << "MolState::Read(): Error. Wrong format in geometry: ["+geom_iStr.str()+"]\n";
      exit(1);
    }
    atoms.push_back(tmp_atom);
  }
  //std::cout << "Read geometry: DONE" << std::endl;
  
  //------------ Convert atomic names to masses --------------------------
  for (i=0; i<NAtoms(); i++) {

    //std::cout << "Atom name=" << getAtom(i).Name().c_str() << std::endl;
    //std::cout << "Value=" << node_amu_table.read_node_double_value(getAtom(i).Name().c_str()) << std::endl;
								   
    getAtom(i).Mass()=node_amu_table.read_node_double_value(getAtom(i).Name().c_str());
  }
  //std::cout << "Read masses: DONE" << std::endl;   
   
  //------------ Read Normal Modes ---------------------------------------
  NormalMode tmp_normMode(NAtoms(),0); // one temp. norm mode

  normModes.clear();

  xml_node node_nmodes(node_state,"normal_modes",0);
  My_istringstream nmodes_iStr(node_nmodes.read_string_value("text"));
  
  for (i=0; i < tmp_nNormMds; i++)
    normModes.push_back( tmp_normMode );
  
  int nModesPerLine=3;  // three number of vib. modes per Line

  int nLines;
  nLines = tmp_nNormMds / nModesPerLine;
  if ( tmp_nNormMds % nModesPerLine != 0)
    nLines++;

  for (k = 0; k < nLines; k++)  // number of blocks with 3 norm.modes. ("lines")
    {
      int current_nModesPerLine = nModesPerLine;
      // for the last entree, nModesPerString may differ from 3
      if (nLines - 1 == k)
  	if ( tmp_nNormMds % nModesPerLine != 0 )
	  current_nModesPerLine = tmp_nNormMds % nModesPerLine;
 
      for (i=0; i < NAtoms(); i++)   
  	for (j=0; j < current_nModesPerLine; j++) 
  	  for (l=0; l < CARTDIM; l++) {
	    normModes[k*nModesPerLine+j].getDisplacement()[i*CARTDIM+l]=nmodes_iStr.getNextDouble();
	    if (nmodes_iStr.fail()) {
	      std::cout<<"MolState::Read(): Error. Wrong format in normal modes: ["+nmodes_iStr.str()+"]\n";
	      exit(1);
	    }
	  }
    }
  
   
  //------------ Read Frequencies ----------------------------------------
  xml_node node_freq(node_state,"frequencies",0);
  My_istringstream freq_iStr(node_freq.read_string_value("text"));
  
  for (i=0; i < tmp_nNormMds; i++)
  {
    getNormMode(i).getFreq()=freq_iStr.getNextDouble();
    if (freq_iStr.fail()) {
      std::cout << "\nError: format error at frequency #"<<i<< " or less than " << tmp_nNormMds <<" frequencies found\n";
      exit(1);
    }
    if (getNormMode(i).getFreq()<=0) {
      std::cout <<"\nError. The frequency ["<<getNormMode(i).getFreq() <<"] is negative\n";
      exit(1);
    }
  }
  //std::cout << "Read frequences: DONE\n" ;
  
  // Now MolState is in a good shape, and some transformations can be performed

  //------------ 1. Un-mass-weight normal modes ----------------------------
  bool if_massweighted=node_nmodes.read_bool_value("if_mass_weighted");


  //------------ 1. mass un-weight normal modes, if needed (QChem-->ACES format;) ----------------
  // qchem if_massweighted="true"; aces if_massweighted="false";
  reduced_masses.Adjust(NNormModes(),1);
  reduced_masses.Set(1);
  if (if_massweighted)
    {
      //std::cout << "Doing mass-weighted coordinates.....\n" ;
      // Read atomic names from "...->normal_modes->atoms"
      My_istringstream atoms_iStr(node_nmodes.read_string_value("atoms"));
      
      //std::cout << "Atoms: " << atoms_iStr.str() << std::endl;
      
      std::vector<Atom> normalModeAtoms;
  
      normalModeAtoms.clear();

      for (i=0; i<NAtoms(); i++) {

	std::string tmp_atomName;
	//Get atomic name:
	atoms_iStr.getNextWord(tmp_atomName);
	tmp_atom.Name() = tmp_atomName;
	normalModeAtoms.push_back(tmp_atom);
      }
      
      // Get masses for each atomic name:
      for (i=0; i<NAtoms(); i++){

	//std::cout << "Atom name=" << normalModeAtoms[i].Name().c_str() << std::endl;
	//std::cout << "Value=" << node_amu_table.read_node_double_value(normalModeAtoms[i].Name().c_str()) << std::endl;
	normalModeAtoms[i].Mass()=node_amu_table.read_node_double_value(normalModeAtoms[i].Name().c_str());
      }
      
      // Mass-un-weight normal modes:
      for (int nm=0; nm<NNormModes(); nm++)
	for (int a=0; a<NAtoms(); a++)
	  for (int i=0; i<CARTDIM; i++ )
	    getNormMode(nm).getDisplacement().Elem1(a*CARTDIM+i) *= sqrt(normalModeAtoms[a].Mass());

      // normalize each normal mode (/sqrt(norm) which is also /sqrt(reduced mass)):
      // KMatrix reduced_masses(NNormModes(),1);
      for (int nm=0; nm<NNormModes(); nm++) {
	reduced_masses[nm] = 0;
	for (int a=0; a<NAtoms(); a++)
	  for (int i=0; i<CARTDIM; i++ )
	    reduced_masses[nm]+= getNormMode(nm).getDisplacement().Elem1(a*CARTDIM+i) *
	      getNormMode(nm).getDisplacement().Elem1(a*CARTDIM+i);
      }
      // reduced_masses.Print("Reduced masses:");
      // Normalize:
      for (int nm=0; nm<NNormModes(); nm++)
	for (int a=0; a<NAtoms(); a++)
	  for (int i=0; i<CARTDIM; i++ )
	    getNormMode(nm).getDisplacement().Elem1(a*CARTDIM+i)/=sqrt(reduced_masses[nm]);

    }
  

  //FIXIT: Separate the code below into a member function of this class.
  
  // ------------ 1.5 Find geometry from the vertial gradient if available ------------
    
  /* Notes on the gradient implementation:
   * Detection of gradient node triggers use of gradient and changes meaning of nodes: geometry, normal modes, frequncies, and excitation_energy. If a gradient node is present in an electronic state, the geometry, normal modes, and frequncies nodes are expected to descibe the initial state. The input excitation energy has to be the vertial excitation energy at the initial state geometry. The target state geometry and adiabatic excitation energy will be calculated using the vertial gradient method.
   * Pawel, May 2022
   * */
  if (node_state.find_subnode("gradient")) {
    std::cout 
      << " State geometry will be calculated with the vertical gradient (VG) approximation." 
      << std::endl;
    xml_node node_gradient(node_state, "gradient", 0);
    // This version supports gradient only in atomic units, a.u., (Q-Chem opt output default)
    std::string units = node_gradient.read_string_value("units");
    if (units != "a.u.") {
      std::cout 
        << "\nError! Gradient reported in unsupported units: \"" 
        << units 
        << "\"\n  (this version of ezFCF supports only \"a.u.\")\n\n";
      exit(1);
    }
    My_istringstream grad_iStr(node_gradient.read_string_value("text"));
    // gradient is stored as a column, i.e, a 3N x 1 matrix 
    gradient.Adjust(CARTDIM * NAtoms(), 1);
    for (int i = 0; i < NAtoms(); i++) {
      // Read the x, y, and z components of each vector
      for (int j = 0; j < CARTDIM; j++) {
        double cartesian_component = grad_iStr.getNextDouble();
        if (grad_iStr.fail()) {
          std::cout 
            << "MolState::Read(): Error. Wrong format in gradient: [" 
            << grad_iStr.str() 
            << "]" 
            << std::endl;
          exit(1);
        }
        gradient.Elem1(CARTDIM * i + j) = cartesian_component;
      }
    }
  }
  /* End of parsing of the vertical gradient node. */
  if (gradient.Size() > 0) {
    // $\Delta = \Omega ^{-2} D ^{-1} M ^{-1/2} g _{(2)} ^{(x)}$

    /* $\Delta$ is a vector of displacement between the target state normal coordinates $q ^{(2)}$ and the initial state normal coordinates $q ^{(1)}$:
     * $$ q ^{(2)} = q ^{(1)} + \Delta $$
     * (in parallel approximation without frequency shifts there is no Duschinsky matrix in the last equation, and as the two sets of normal modes are parallel, $\Delta$ is also the shift between the two geometries in the normal mode coordinates)
     * $\Omega$ is a diagonal matrix (3N - 5/6 x 3N - 5/6) of harmonic frequencies 
     * $D$ has normal modes as its columns. $D$ is a rectangular matrix (3N x 3N - 5/6) diagonalizing the mass-weighted Hessian matrix:
     * $$ H = D \Omega ^2 D ^{-1} $$
     * One dimension of $D$ is shorter because the normal coordinates coresponding to zero frequency modes do not contribute to vibrational spectrum but descibe translations and rotations. 
     * M is a diagonal matrix of dim 3N x 3N of atomic masses
     * g _(2) ^{(x)} is the gradient of the target state energy surface calculated at the geometry of the initial state in cartesian (non-mass-weighted) coordinates (${}^{(x)}$). The number ${}_(2)$ indicates that the value was calculated at the second (target) state.
     * 
     * Additionally, to find $R _e ^(2)$ (the equilibrium geometry of the second (target) state) from $\Delta$
     * $$ R _e ^(2) = R _e ^{(1)} - M ^{-1/2} D \Delta $$
     * or in terms of a gradient
     * $$ R _e ^(2) = R _e ^{(1)} - M ^{-1/2} D \Omega ^{-2} D ^{-1} M ^{-1/2} g _{(2)} ^{(x)} $$
     *
     * $M$ is in the $-1/2$ power in both cases as the right most is a transformation of the **gradient** from cartesian to mass-weighted coordinates, while the leftmost is a transfromation of **coordinates** from a mass-weighted vector to a mass-un-weighed vector.
     *
     * $\Delta$ as well as the cartesian displacement ($M ^{-1/2} D \Delta$) are expressed in a.u. and are converted to Angstroms only just before addition to the geometry. For this reason both matrices ($\Omega$ and $M$) and the gradient vector have values corresponding to the the units of a.u.
     */

    //FIXIT: This needs to be cleaned up a bit in the future: use AU consistently,
    //Move some pieces into functions of appropriate classes (e.g., normalmode)
    KMatrix mass_matrix_minus_half(CARTDIM * NAtoms(), CARTDIM * NAtoms(), true);
    for (int i = 0; i < NAtoms(); i++) {
      double mass = atoms[i].Mass() * AMU_2_ELECTRONMASS;
      double mass_inv = 1.0 / mass;
      double sqrt_mass_inv = sqrt(mass_inv);
      for (int j = 0; j < CARTDIM; j++) {
        mass_matrix_minus_half.Elem2(CARTDIM * i + j, CARTDIM * i + j) = sqrt_mass_inv;
      }
    }

    KMatrix d_matrix(CARTDIM * NAtoms(), NNormModes(), true);
    KMatrix Omega_matrix_minus2(NNormModes(), NNormModes(), true);
    for (int j = 0; j < NNormModes(); j++){
      double freq = normModes[j].getFreq();
      freq *= WAVENUMBERS2EV; // first to eV
      freq *= EV2HARTREE; // then eV to a.u.
      Omega_matrix_minus2.Elem2(j, j) = 1 / freq / freq;
      for (int i = 0; i < CARTDIM * NAtoms(); i++){
        d_matrix.Elem2(i, j) = normModes[j].getDisplacement()[i];
      }
    }
    
    KMatrix delta(gradient, true); // true initilizes values of delta with values of gradient
    delta.LeftMult(mass_matrix_minus_half);
    delta.LeftMult(d_matrix, true); // true transposes the matrix, D ^T = D ^{-1}
    delta.LeftMult(Omega_matrix_minus2);

    // R _e ^{(2)} = R _e ^{(1)} - M ^{-1/2} D \Delta
    KMatrix geometry_shift(delta, true); // true for copying the data
    geometry_shift.LeftMult(d_matrix);
    geometry_shift.LeftMult(mass_matrix_minus_half);

    geometry_shift *= AU2ANGSTROM;

    for (int i = 0; i < NAtoms(); i++) {
      for (int j = 0; j < CARTDIM; j++){
        atoms[i].Coord(j) -= geometry_shift.Elem2(CARTDIM * i + j, 0);
      }
    }
    std::cout << "Target-state geometry calculated with vertical gradient apprx-n:" << std::endl;
    printGeometry();

    for (int i = 0; i < NNormModes(); i++) {
      energy -= 0.5 * delta.Elem1(i) * delta.Elem1(i) / Omega_matrix_minus2.Elem2(i,i) / EV2HARTREE; 
    }

    std::cout 
      << "Adiabatic excitation energy (within VG) = "
      << energy << " eV " << std::endl;
  }
  //end of computing geometry and adiabatic energy by using VG

  //------------ 2. Align geometry if requested ----------------------------
  if_aligned_manually=false;
  double man_rot_x, man_rot_y, man_rot_z;
  Vector3D man_shift;

  size_t manual_coord_transform=node_state.find_subnode("manual_coordinates_transformation");
  
  for (int i = 0; i < manual_coord_transform; ++i) {

    //std::cout  << "Do manual transformation" << std::endl;
    xml_node manual_coord_transform(node_state,"manual_coordinates_transformation", i);
    
    man_rot_x=manual_coord_transform.read_double_value("rotate_around_x"); 
    man_rot_y=manual_coord_transform.read_double_value("rotate_around_y"); 
    man_rot_z=manual_coord_transform.read_double_value("rotate_around_z"); 

    man_shift.getCoord(0)=-manual_coord_transform.read_double_value("shift_along_x"); 
    man_shift.getCoord(1)=-manual_coord_transform.read_double_value("shift_along_y"); 
    man_shift.getCoord(2)=-manual_coord_transform.read_double_value("shift_along_z"); 
       
    std::cout << "Molecular structure and normal modes of this electronic state\nwill be transformed as requested in the input.\n";

    shiftCoordinates(man_shift);
    rotate(man_rot_x*PI, man_rot_y*PI, man_rot_z*PI);
    applyCoordinateThreshold(COORDINATE_THRESHOLD);

    std::cout << "Molecule was shifted by " 
      << man_shift[0] << ", "
      << man_shift[1] << ", "
      << man_shift[2] << " in x, y, and z." 
      << std::endl;

    std::cout << "Also rotated by " 
      << man_rot_x <<"*pi, " 
      << man_rot_y <<"*pi, and " 
      << man_rot_z <<"*pi around x, y, and z axes.\n";

    // Printing difference from ground state is not implemented:
    // This function would require knowldege about the ground state, and 
    // as this is a Read state function it appears inappropriet to add it here.
    // It appear approprieate, however, to move this manual rotation out of the 
    // Read function (same as the vertical gradient) and apply the post processing 
    // once all input states are parsed.
    /* double diff_from_ground = this->getGeomDifference(*this); */
    /* std::cout << "The norm of the geometry difference from the initial state is " << sqrt(diff_from_ground) <<"\n"; */

    if_aligned_manually=true;
  }

  
  //------------ 3. Reorder normal modes if requested --------------------

  size_t reorder_nmodes=node_state.find_subnode("manual_normal_modes_reordering");
  if ( reorder_nmodes ) {
    
    xml_node node_nmodes_reorder(node_state,"manual_normal_modes_reordering",0);
    My_istringstream reorder_iStr(node_nmodes_reorder.read_string_value("new_order"));
      
    if_nm_reordered_manually=false;//?? 
    normModesOrder.clear();

    std::cout << "New normal modes order was requested:\n" << reorder_iStr.str() <<"\n";
    int tmpInt;
    for (int nm=0; nm < NNormModes(); nm++)
      {
	tmpInt=reorder_iStr.getNextInt();
	
	//input error check:
	if (reorder_iStr.fail()) {
	  std::cout << "\nFormat error: non-numeric symbol or less entries then the number of normal modes\n\n";
	}
	if ( (tmpInt<0) or (tmpInt>=NNormModes()) ) {
	  std::cout << "\nError: normal mode number ["<< tmpInt<<"] is out of range [0.."<<NNormModes()-1<<"].\n\n";
	}
	normModesOrder.push_back(tmpInt);
      }
    
    // check if there are duplicates in the list:
    std::vector<int> tmpIntVector, tmpIntVector2;
    tmpIntVector = normModesOrder;
    std::sort( tmpIntVector.begin(), tmpIntVector.end() );
    tmpIntVector2 = tmpIntVector;
    std::vector<int>::const_iterator intVec_iter;
    intVec_iter= unique( tmpIntVector.begin(), tmpIntVector.end() );
    if (intVec_iter != tmpIntVector.end())
    {
      std::cout << "\nFormat error: there are non unique entries. Check the sorted list:\n";
      for (std::vector<int>::const_iterator tmp_iter=tmpIntVector2.begin(); tmp_iter!=tmpIntVector2.end(); tmp_iter++)
        std::cout << ' ' << *tmp_iter;
      std::cout<<'\n';
    }
      
    // backup normal modes
    std::vector<NormalMode> oldNormModes;
    oldNormModes.clear();
    for (int nm=0; nm < NNormModes(); nm++) {
      tmp_normMode = getNormMode(nm);
      oldNormModes.push_back( tmp_normMode );
    }
    
    // copy normal modes using new order
    for (int nm=0; nm < NNormModes(); nm++)
      getNormMode(nm) = oldNormModes[  normModesOrder[nm] ];
    
    std::cout << "Normal modes were reordered accordingly.\n";
    if_nm_reordered_manually=true;
    
  }
  else
    for (int nm=0; nm < NNormModes(); nm++)
      normModesOrder.push_back(nm);


  //------------ 4. Reorder atoms if requested --------------------

  Atom tmp_Atom; // one temp. norm mode
  
  size_t reorder_atoms=node_state.find_subnode("manual_atoms_reordering");
   
  if (reorder_atoms) {

    xml_node reorder_atoms(node_state,"manual_atoms_reordering",0);
    My_istringstream reorder_iStr(reorder_atoms.read_string_value("new_order"));
    std::cout << "New order of atoms was requested:\n" << reorder_iStr.str() <<"\n";
    
    std::vector<int> atomsOrder;
    atomsOrder.clear();
    
    int tmpInt;
    for (int nm=0; nm < NAtoms(); nm++) {

      tmpInt=reorder_iStr.getNextInt();
      //input error check:
      if (reorder_iStr.fail()) {
	
	std::cout << "\nFormat error: non numeric symbol or less entries then the number of atoms\n\n";
	exit(1);
      }
      if ( (tmpInt<0) or (tmpInt>=NAtoms()) ) {
	std::cout << "\nError: stom number ["<< tmpInt<<"] is out of range [0.."<<NAtoms()-1<<"].\n";
	exit(1);
      }
      atomsOrder.push_back(tmpInt);
    }
      
    // check if there are duplicates in the list:
    std::vector<int> tmpIntVector, tmpIntVector2;
    tmpIntVector = atomsOrder;
    std::sort( tmpIntVector.begin(), tmpIntVector.end() );
    tmpIntVector2 = tmpIntVector;
    std::vector<int>::const_iterator intVec_iter;
    intVec_iter= unique( tmpIntVector.begin(), tmpIntVector.end() );
    if (intVec_iter != tmpIntVector.end()) {
      std::cout << "\nFormat error: there are non unique entries. Check the sorted list:\n";
      for (std::vector<int>::const_iterator tmp_iter=tmpIntVector2.begin(); tmp_iter!=tmpIntVector2.end(); tmp_iter++)
	std::cout << ' ' << *tmp_iter << std::endl;

      exit(1);
    }
    
    // backup molecular geometry and normal modes
    std::vector<Atom> oldAtoms;
    oldAtoms.clear();
    for (int a=0; a < NAtoms(); a++) {

      tmp_Atom = getAtom(a);
      oldAtoms.push_back( tmp_Atom );
    }
    std::vector<NormalMode> oldNormModes;
    oldNormModes.clear();
    for (int nm=0; nm < NNormModes(); nm++) {

      tmp_normMode = getNormMode(nm);
      oldNormModes.push_back( tmp_normMode );
    }
      
    // copy molecular geometry using the new order of atoms
    for (int a=0; a < NAtoms(); a++)
      getAtom(a) = oldAtoms[  atomsOrder[a] ];
    
    // copy normal modes using the new order of atoms
    for (int nm=0; nm < NNormModes(); nm++)
      for (int a=0; a < NAtoms(); a++) 
	for (int k=0; k < CARTDIM; k++)
	  getNormMode(nm).getDisplacement()[a*CARTDIM+k]=oldNormModes[nm].getDisplacement()[  atomsOrder[a]*CARTDIM + k  ];
    
    std::cout << "Atoms were reordered accordingly.\n";
  }

  return true;
}
