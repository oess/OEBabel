/**********************************************************************
Copyright (C) 2005, 2006, 2007, 2008, 2009, 2010, 2013 OpenEye Scientific Software, Inc.
***********************************************************************/
#include "openeye.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "oeplatform.h"
#include "oesystem.h"
#include "oedepict.h"

#include "babel.itf"

using namespace OEPlatform;
using namespace OESystem;
using namespace OEChem;
using namespace std;

static const string BABEL_VERSION = "3.4beta";
static const string PROGNAME = "babel";

// Tracks the number of molecule in and out
static unsigned int ICOUNT = 0;

template <class MolType>
class OENonConstMolFunc
{
public:
  virtual ~OENonConstMolFunc() {};
  virtual bool operator()(MolType &mol) = 0;
};

////////////////////
// Hydrogen handling
////////////////////
template<class MolType>
class OEHydrogens : public OENonConstMolFunc<MolType>
{
  unsigned int _hydrogens;
  enum { Undefined, Add, AddMDLImplicit, Delete };
public:
  OEHydrogens(const string hydrogens)
    : OENonConstMolFunc<MolType>(), _hydrogens(Undefined)
  {
    if (hydrogens=="add")
      _hydrogens = Add;
    else if (hydrogens=="add_mdl_implicit")
      _hydrogens = AddMDLImplicit;
    else if (hydrogens=="delete")
      _hydrogens = Delete;
    else
      OEThrow.Error("Illegal hydrogen action requested.");
  }

  bool operator()(MolType &mol)
  {
    bool ret = false;
    switch(_hydrogens)
    {
      case Add:
        ret = OEAddExplicitHydrogens(mol);
        if (OEHasResidues(mol))
          OEResidueHydrogens(mol);
        return ret;
      case AddMDLImplicit:
        OEAssignMDLHydrogens(mol);
        return true;
      case Delete:
        return OESuppressHydrogens(mol);
      default:
        OEThrow.Error("Illegal hydrogen action requested.");
        return false;
    }
  }
};

/////////////////
// Invalid stereo
/////////////////
class OEDeleteInvalidStereo : public OENonConstMolFunc<OEMolBase>
{
public:
  bool operator ()(OEMolBase &mol)
  {
    vector<OEAtomBase *> nbrs;
    unsigned int ia=0,ib=0;

    OEPerceiveChiral(mol);
    for (OEIter<OEAtomBase> atom=mol.GetAtoms();atom;++atom)
    {
      if (atom->HasStereoSpecified() && !(atom->IsChiral()))
      {
        nbrs.clear();
        for (OEIter<OEAtomBase> nbr=atom->GetAtoms();nbr;++nbr)
          nbrs.push_back(atom);
        atom->SetStereo(nbrs,OEAtomStereo::Tetrahedral,OEAtomStereo::Undefined);
        ++ia;
      }
    }
    for (OEIter<OEBondBase> bond=mol.GetBonds();bond;++bond)
    {
      if (bond->GetOrder()!=2)
        continue;
      if (bond->HasStereoSpecified() && !(bond->IsChiral()))
      {
        nbrs.clear();
        nbrs.push_back(bond->GetBgn());
        nbrs.push_back(bond->GetEnd());
        bond->SetStereo(nbrs,OEBondStereo::CisTrans,OEBondStereo::Undefined);
        ++ib;
      }
    }
    if (ia > 0)
      OEThrow.Warning("invalid stereo deleted for %d atoms: %s",ia,mol.GetTitle());
    if (ib > 0)
      OEThrow.Warning("invalid stereo deleted for %d bonds: %s",ib,mol.GetTitle());
    return true;
  }
};

/////////////////
// Stereo from 3D
/////////////////
static bool IsMDLConnectionTable(unsigned int fmt)
{
  switch (fmt)
  {
    case OEFormat::MDL:
    case OEFormat::SDF:
      return true;
    default:
      return false;
  }
}

class OEStereoFrom3D : public OENonConstMolFunc<OEMolBase>
{
  unsigned int _ofmt;
public:
  OEStereoFrom3D(unsigned int ofmt) : _ofmt(ofmt) {}

  bool operator ()(OEMolBase &mol)
  {
    bool ret = false;
    ret = mol.SetDimension(3);
    ret = OE3DToInternalStereo(mol) && ret;
    if (IsMDLConnectionTable(_ofmt))
      OEMDLPerceiveParity(mol);
    return ret;
  }
};

//////////////////////////
// MDL correct bond stereo
//////////////////////////
class OEMDLCorrectStereo : public OENonConstMolFunc<OEMolBase>
{
  unsigned int _ifmt;
public:
  OEMDLCorrectStereo(unsigned int ifmt) : _ifmt(ifmt)
  {
    if (!IsMDLConnectionTable(_ifmt))
      OEThrow.Warning("Will only correct MDL stereo when reading from MDL files");
  }

  bool operator ()(OEMolBase &mol)
  {
    if (IsMDLConnectionTable(_ifmt) &&
        mol.GetDimension()==2       &&
        OEMDLHasIncorrectBondStereo(mol))
    {
      bool ret = false;
      ret = OEMDLCorrectBondStereo(mol);
      ret = OEMDLStereoFromBondStereo(mol) && ret;
      return ret;
    }
    return true;
  }
};

/////////////
// OESD2Title
/////////////
class OESD2Title : public OENonConstMolFunc<OEMolBase>
{
  string _tag;
public:
  OESD2Title(string tag) : _tag(tag) {}

  bool operator ()(OEMolBase &mol)
  {
    if (OEHasSDData(mol, _tag))
      return mol.SetTitle(OEGetSDData(mol, _tag));

    OEThrow.Warning("%s missing SDTag %s", mol.GetTitle(), _tag.c_str());
    return false;
  }
};

///////////////
// OESD2TitleMC
///////////////
class OESD2TitleMC : public OENonConstMolFunc<OEMCMolBase>
{
  string _tag;
  OESD2Title _scmolfunc;
public:
  OESD2TitleMC(string tag) : _tag(tag), _scmolfunc(tag) {}

  bool operator ()(OEMCMolBase &mol)
  {
    bool molHasTitle = OEHasSDData(mol,_tag);
    string molTitle;
    bool ret = true;
    if (molHasTitle)
    {
      molTitle = OEGetSDData(mol,_tag);
      ret = _scmolfunc((OEMolBase&)mol) && ret;
    }

    for (OEIter<OEConfBase> conf=mol.GetConfs();conf;++conf)
    {
      if (OEHasSDData(*conf,_tag))
        ret = conf->SetTitle(OEGetSDData(*conf,_tag)) && ret;
      else if (molHasTitle)
        ret = conf->SetTitle(molTitle) && ret;
      else
      {
        OEThrow.Warning("%s missing SDTag %s", mol.GetTitle(), _tag.c_str());
        return false;
      }
    }
    return ret;
  }
};

////////////////////
// Perceive residues
////////////////////
class OEPerceiveResiduesFunc : public OENonConstMolFunc<OEMolBase>
{
  unsigned int _mask;
public:
  OEPerceiveResiduesFunc(unsigned int mask) : _mask(mask) {}

  bool operator ()(OEMolBase &mol)
  {
    OEPerceiveResidues(mol,_mask);
    return true;
  }
};

///////////////////
// Canonical Kekule
///////////////////
class OECanonicalKekule : public OENonConstMolFunc<OEMolBase>
{
public:
  bool operator ()(OEMolBase &mol)
  {
    OEFindRingAtomsAndBonds(mol);
    OEAssignAromaticFlags(mol, OEAroModelOpenEye);

    bool ret = true;
    for (OEIter<OEBondBase> bond=mol.GetBonds(OEIsAromaticBond());bond;++bond)
      ret = bond->SetIntType(5) && ret;

    OECanonicalOrderAtoms(mol);
    OECanonicalOrderBonds(mol);
    OEClearAromaticFlags(mol);
    ret = OEKekulize(mol) && ret;
    return ret;
  }
};

////////////////////////
// 2D depict coordinates
////////////////////////

class DepictCoordinates : public OENonConstMolFunc<OEMolBase>
{
  unsigned int _ofmt;
public:
  DepictCoordinates(unsigned int ofmt) : _ofmt(ofmt) {}

  bool operator ()(OEMolBase &mol)
  {
    bool ret = false;
    ret = OEDepict::OEAddDepictionHydrogens(mol) && ret;
    ret = OEDepict::OEDepictCoordinates(mol) && ret;
    if (IsMDLConnectionTable(_ofmt))
        ret = OEMDLPerceiveBondStereo(mol) && ret;
    return mol.SetDimension(2) && ret;
  }
};

////////////////////////
// Output titles to a file
////////////////////////
class OEOutputNames : public OENonConstMolFunc<OEMolBase>
{
  oeofstream _ofs;
public:
  OEOutputNames(string fname)
  {
    if (!_ofs.open(fname))
      OEThrow.Fatal("Unable to open namefile: " + fname);
  }

  bool operator ()(OEMolBase &mol)
  {
    _ofs << mol.GetTitle() << oeendl;
    return true;
  }
};

/////////////////////////
// Count input conformers
/////////////////////////
class OEConfCountBase
{
protected:
  unsigned int _count;
public:
  OEConfCountBase() : _count(0) {}
  unsigned int GetCount() const { return _count; }
};

class OEConfCount : public OEConfCountBase,
                    public OENonConstMolFunc<OEMolBase>
{
public:
  bool operator ()(OEMolBase &)
  {
    _count++;
    return true;
  }
};

////////////////////////////////////
// Count input conformers for MCMols
////////////////////////////////////
class OEConfCountMC : public OEConfCountBase,
                      public OENonConstMolFunc<OEMCMolBase>
{
public:
  bool operator ()(OEMCMolBase &mol)
  {
    _count += mol.NumConfs();
    return true;
  }
};

//////////////////////////////////////////
// Abstraction for the molecule output end
//////////////////////////////////////////
class OECountOutputBase
{
protected:
  unsigned int _ocount;
  unsigned int _oconfcount;
public:
  OECountOutputBase() : _ocount(0), _oconfcount(0) {}
  virtual ~OECountOutputBase() {};

  unsigned int GetCount() const { return _ocount; }

  void AddConfCount(OEMolBase &) { _oconfcount++; }
  void AddConfCount(OEMCMolBase &mol) { _oconfcount+=mol.NumConfs(); }
  unsigned int GetConfCount() const { return _oconfcount; }
};

template <class StreamType, class MolType>
class OEOutputMolFunc : public OECountOutputBase
{
public:
  virtual unsigned int operator()(StreamType &ofs, MolType &mol) = 0;
};

////////////////////////////////////////////////////
// Dump the molecule out normally
////////////////////////////////////////////////////
template <class StreamType, class MolType>
class OENormalWriteMolecule : public OEOutputMolFunc<StreamType, MolType>
{
  using OEOutputMolFunc<StreamType, MolType>::_ocount;
public:
  unsigned int operator()(StreamType &ofs, MolType &mol)
  {
    OEWriteMolecule(ofs, mol);
    _ocount++;
    this->AddConfCount(mol);
    return 1;
  }
};

template<class T>
class MolProxy
{
};

template<>
class MolProxy<OEMCMolBase> : public OEMol
{
};

template<>
class MolProxy<OEMolBase> : public OEGraphMol
{
};

////////////////////////////////////////////////////
// Split the molecule into it's connected components
////////////////////////////////////////////////////
template <class StreamType, class MolBaseT>
class OEParts2Mols : public OEOutputMolFunc<StreamType, MolBaseT>
{
  using OEOutputMolFunc<StreamType, MolBaseT>::_ocount;
public:
  unsigned int operator()(StreamType &ofs, MolBaseT &mol)
  {
    unsigned int maxidx = mol.GetMaxAtomIdx();
    OEMallocaPtr<unsigned int> parts = OEMalloca(sizeof(unsigned int)*maxidx);
    unsigned int count = OEDetermineComponents(mol,parts);

    MolProxy<MolBaseT> partmol;
    OEPartPred pred(parts, maxidx);
    for (unsigned int i=1; i<=count; i++)
    {
      pred.SelectPart(i);
      OESubsetMol(partmol,mol,pred);
      OEWriteMolecule(ofs,(MolBaseT &)partmol);
      this->AddConfCount((MolBaseT &)partmol);
    }
    _ocount+=count;

    return count;
  }
};

/////////////////////////////////////
// Single-conformerize a MCMol stream
////////////////////////////////////.
template <class StreamType>
class OEMC2SC : public OEOutputMolFunc<StreamType, OEMCMolBase>
{
  using OEOutputMolFunc<StreamType, OEMCMolBase>::_ocount;
  using OEOutputMolFunc<StreamType, OEMCMolBase>::AddConfCount;
public:
  unsigned int operator()(StreamType &ofs, OEMCMolBase &mol)
  {
    OEWriteMolecule(ofs, *mol.GetActive());
    _ocount++;
    AddConfCount(*mol.GetActive());
    return 1;
  }
};

namespace OutputFuncTypes
{
const unsigned int Normal     = 0;
const unsigned int Parts2Mols = 1;
const unsigned int MC2SC      = 2;
}

//////////////////////////////////
// Container for the mol functions
//////////////////////////////////
class OEBabel
{
  typedef vector<OENonConstMolFunc<OEMolBase>* > SCContainer;
  typedef SCContainer::iterator SCIterator;
  typedef vector<OENonConstMolFunc<OEMCMolBase>* > MCContainer;
  typedef MCContainer::iterator MCIterator;
  SCContainer _scfuncs;
  MCContainer _mcfuncs;
  OECountOutputBase *_output;
  unsigned int _outtype;
  OEHeader _header;
  OEBabel();
public:
  OEBabel(int argc, char *argv[], const OEInterface &itf)
    : _output(0),
      _outtype(0),
      _header(argc, argv, PROGNAME.c_str(), BABEL_VERSION.c_str(), "", itf)
  {}
  ~OEBabel()
  {
    for (SCIterator func=_scfuncs.begin(); func!=_scfuncs.end(); ++func)
      delete *func;

    for (MCIterator func=_mcfuncs.begin(); func!=_mcfuncs.end(); ++func)
      delete *func;
  }

  OEHeader &GetHeader()
  {
    return _header;
  }

  void AddTransform(OENonConstMolFunc<OEMolBase> *func)
  {
    _scfuncs.push_back(func);
  }

  void AddTransform(OENonConstMolFunc<OEMCMolBase> *func)
  {
    _mcfuncs.push_back(func);
  }

  void ApplyTransform(OEMolBase &mol)
  {
    for (SCIterator func=_scfuncs.begin(); func!=_scfuncs.end(); ++func)
      if (!(**func)(mol))
        OEThrow.Warning("Problem processing %s", mol.GetTitle());
  }

  void ApplyTransform(OEMCMolBase &mol)
  {
    for (MCIterator func=_mcfuncs.begin(); func!=_mcfuncs.end(); ++func)
      if (!(**func)(mol))
        OEThrow.Warning("Problem processing %s", mol.GetTitle());

    // Call the single-conf functions as well
    ApplyTransform(static_cast<OEMolBase &>(mol));
  }


  void SetOutputType(unsigned int t) { _outtype = t; }
  unsigned int GetOutputType() { return _outtype; }

  template<class StreamType>
  OEOutputMolFunc<StreamType,OEMolBase> *GetOutputFunc(OEMolBase &)
  {
    switch(GetOutputType())
    {
      case OutputFuncTypes::Normal:
        return new OENormalWriteMolecule<StreamType, OEMolBase>;
      case OutputFuncTypes::Parts2Mols:
        return new OEParts2Mols<StreamType, OEMolBase>;
      default:
        OEThrow.Fatal("Unknown output type %u", GetOutputType());
        return 0;
    }
  }
  template<class StreamType>
  OEOutputMolFunc<StreamType,OEMCMolBase> *GetOutputFunc(OEMCMolBase &)
  {
    switch(GetOutputType())
    {
      case OutputFuncTypes::Normal:
        return new OENormalWriteMolecule<StreamType, OEMCMolBase>;
      case OutputFuncTypes::Parts2Mols:
        return new OEParts2Mols<StreamType, OEMCMolBase>;
      case OutputFuncTypes::MC2SC:
        return new OEMC2SC<StreamType>;
      default:
        OEThrow.Fatal("Unknown output type %u", GetOutputType());
        return 0;
    }
  }

  template <class StreamType, class MolBaseT>
  unsigned int WriteMolecule(StreamType &ofs, MolBaseT &mol)
  {
    static OEOutputMolFunc<StreamType, MolBaseT> *ofunc =
      GetOutputFunc<StreamType>(mol);
    if (!_output)
      _output = ofunc;

    if (_output != ofunc)
      OEThrow.Fatal("Two different output methods are being used!");

    return (*ofunc)(ofs, mol);
  }

  unsigned int GetOCount() const
  {
    if (_output)
      return _output->GetCount();
    return 0;
  }
  unsigned int GetOConfCount() const
  {
    if (_output)
      return _output->GetConfCount();
    return 0;
  }
};

/////////////////////////////////////////////////////////////////////////////
static bool IsMCMolInput(const OEInterface &itf, unsigned int ifmt)
{
  if ((itf.Get<bool>("-mc") ||
       itf.Get<bool>("-mc_isomer") ||
       ifmt==OEFormat::OEB) &&
      !itf.Get<bool>("-sc"))
    return true;
  return false;
}

/////////////////////////////////////////////////////////////////////////////

static bool Is3DFormat(unsigned int oefmt)
{
  switch (oefmt)
  {
    case (OEFormat::OEB):
    case (OEFormat::SDF):
    case (OEFormat::MDL):
    case (OEFormat::MOL2):
    case (OEFormat::MOL2H):
    case (OEFormat::XYZ):
    case (OEFormat::PDB):
    case (OEFormat::MOPAC):
    case (OEFormat::MMOD):
    case (OEFormat::RDF):
      return true;
    default:
      return false;
  }
}

/////////////////////////////////////////////////////////////////////////////
static unsigned int PerceiveResiduesMask(const OEInterface &itf)
{
  unsigned int mask = 0;
  if (itf.Get<bool>("-perceive_residues_preserve_ChainID"))
    mask |= OEPreserveResInfo::ChainID;
  if (itf.Get<bool>("-perceive_residues_preserve_ResidueNumber"))
    mask |= OEPreserveResInfo::ResidueNumber;
  if (itf.Get<bool>("-perceive_residues_preserve_ResidueName"))
    mask |= OEPreserveResInfo::ResidueName;
  if (itf.Get<bool>("-perceive_residues_preserve_SerialNumber"))
    mask |= OEPreserveResInfo::SerialNumber;
  if (itf.Get<bool>("-perceive_residues_preserve_AlternateLocation"))
    mask |= OEPreserveResInfo::AlternateLocation;
  if (itf.Get<bool>("-perceive_residues_preserve_InsertCode"))
    mask |= OEPreserveResInfo::InsertCode;
  if (itf.Get<bool>("-perceive_residues_preserve_HetAtom"))
    mask |= OEPreserveResInfo::HetAtom;
  if (itf.Get<bool>("-perceive_residues_preserve_All"))
    mask |= OEPreserveResInfo::All;
  return mask;
}

///////////////////
// Helper functions
string GetOutputFile(OEInterface &itf)
{
  if (itf.Has<string>("-ofmt"))
  {
    unsigned int ofmt=OEGetFileType(itf.Get<string>("-ofmt").c_str());
    if (ofmt==OEFormat::UNDEFINED)
      OEThrow.Fatal("Cannot parse -ofmt: %s",itf.Get<string>("-ofmt").c_str());

    string ofname = itf.Get<string>("-ofmt");
    if (itf.Get<bool>("-ogz"))
      ofname += ".gz";

    return ofname;
  }
  else if (itf.Has<string>("-out"))
    return itf.Get<string>("-out");

  return "";
}

unsigned int GetOutputFormat(OEInterface &itf)
{
  return OEGetFileType(OEGetFileExtension(GetOutputFile(itf).c_str()));
}

bool GetBool(OEInterface &itf, const char *name)
{
  return itf.Get<bool>(name);
}

bool HasInt(OEInterface &itf, const char *name)
{
  return itf.Has<int>(name);
}

int GetInt(OEInterface &itf, const char *name)
{
  return itf.Get<int>(name);
}

bool HasString(OEInterface &itf, const char *name)
{
  return itf.Has<string>(name);
}

string GetString(OEInterface &itf, const char *name)
{
  return itf.Get<string>(name);
}

//////////////////////////////////////////////////////////
/// OpenFormattedOStream - open and format the input stream
//////////////////////////////////////////////////////////
template<class StreamType>
bool OpenFormattedOStream(StreamType &ofs, OEInterface &itf,
                          string &ofname, OEBabel *babel=0)
{
  unsigned int mdloflags=0;
  unsigned int mfoflags=0;
  unsigned int mol2oflags=0;
  unsigned int mopacoflags=0;
  unsigned int pdboflags=0;
  unsigned int smioflags=0;
  unsigned int mmodoflags=0;
  unsigned int oflags=0;

  if (!ofs.open(ofname))
    OEThrow.Fatal("Unabled to open %s for output", ofname.c_str());

  unsigned int ofmt = ofs.GetFormat();

  if (babel && ofmt == OEFormat::OEB)
    OEWriteHeader(ofs, babel->GetHeader());

  if (GetBool(itf, "-canonical_isomeric_smiles"))
  {
    if (ofmt!=OEFormat::SMI && ofmt!=OEFormat::ISM && ofmt!=OEFormat::CAN && ofmt!=OEFormat::USM)
      OEThrow.Fatal("-canonical_isomeric_smiles requires .smi, .can .ism or .usm output.");
    ofs.SetFormat(OEFormat::ISM);
  }

  if (GetBool(itf, "-mdloMCHG"))      mdloflags|=OEOFlavor::MDL::MCHG;
  if (GetBool(itf, "-mdloMISO"))      mdloflags|=OEOFlavor::MDL::MISO;
  if (GetBool(itf, "-mdloMRGP"))      mdloflags|=OEOFlavor::MDL::MRGP;
  if (GetBool(itf, "-mdloMV30"))      mdloflags|=OEOFlavor::MDL::MV30;
  if (GetBool(itf, "-mdloMDLParity")) mdloflags|=OEOFlavor::MDL::MDLParity;
  if (GetBool(itf, "-mdloNoParity"))  mdloflags|=OEOFlavor::MDL::NoParity;
  if (GetBool(itf, "-mdloCurrentParity")) mdloflags|=OEOFlavor::MDL::CurrentParity;
  if (GetBool(itf, "-mdloMMask"))     mdloflags|=OEOFlavor::MDL::MMask;
  if (GetBool(itf, "-mdloPMask"))     mdloflags|=OEOFlavor::MDL::PMask;
  if (GetBool(itf, "-mdloDefault")) mdloflags|=OEOFlavor::MDL::Default;
  if (GetBool(itf, "-pdboBONDS"))   pdboflags|=OEOFlavor::PDB::BONDS;
  if (GetBool(itf, "-pdboORDERS"))  pdboflags|=OEOFlavor::PDB::ORDERS;
  if (GetBool(itf, "-pdboBOTH"))    pdboflags|=OEOFlavor::PDB::BOTH;
  if (GetBool(itf, "-pdboCHARGE"))  pdboflags|=OEOFlavor::PDB::CHARGE;
  if (GetBool(itf, "-pdboRADIUS"))  pdboflags|=OEOFlavor::PDB::RADIUS;
  if (GetBool(itf, "-pdboDELPHI"))  pdboflags|=OEOFlavor::PDB::DELPHI;
  if (GetBool(itf, "-pdboTER"))     pdboflags|=OEOFlavor::PDB::TER;
  if (GetBool(itf, "-pdboOEResidues"))     pdboflags|=OEOFlavor::PDB::OEResidues;
  if (GetBool(itf, "-pdboNoResidues"))     pdboflags|=OEOFlavor::PDB::NoResidues;
  if (GetBool(itf, "-pdboCurrentResidues")) pdboflags|=OEOFlavor::PDB::CurrentResidues;
  if (GetBool(itf, "-pdboOrderAtoms")) pdboflags|=OEOFlavor::PDB::OrderAtoms;
  if (GetBool(itf, "-pdboELEMENT"))    pdboflags|=OEOFlavor::PDB::ELEMENT;
  if (GetBool(itf, "-pdboFormalCrg"))  pdboflags|=OEOFlavor::PDB::FormalCrg;
  if (GetBool(itf, "-pdboHETBONDS"))   pdboflags|=OEOFlavor::PDB::HETBONDS;
  if (GetBool(itf, "-pdboDefault")) pdboflags|=OEOFlavor::PDB::Default;
  if (GetBool(itf, "-mopacoXYZ"))     mopacoflags|=OEOFlavor::MOPAC::XYZ;
  if (GetBool(itf, "-mopacoCHARGES")) mopacoflags|=OEOFlavor::MOPAC::CHARGES;
  if (GetBool(itf, "-mopacoDefault")) mopacoflags|=OEOFlavor::MOPAC::Default;
  if (GetBool(itf, "-mfoTitle")) mfoflags|=OEOFlavor::MF::Title;
  if (GetBool(itf, "-mfoDefault")) mfoflags|=OEOFlavor::MF::Default;
  if (GetBool(itf, "-mmodoAtomTypes")) mmodoflags|=OEOFlavor::MMOD::AtomTypes;
  if (GetBool(itf, "-mol2oAtomTypeNames")) mol2oflags|=OEOFlavor::MOL2::AtomTypeNames;
  if (GetBool(itf, "-mol2oBondTypeNames")) mol2oflags|=OEOFlavor::MOL2::BondTypeNames;
  if (GetBool(itf, "-mol2oAtomNames")) mol2oflags|=OEOFlavor::MOL2::AtomNames;
  if (GetBool(itf, "-mol2oOrderAtoms")) mol2oflags|=OEOFlavor::MOL2::OrderAtoms;
  if (GetBool(itf, "-mol2oHydrogens")) mol2oflags|=OEOFlavor::MOL2::Hydrogens;
  if (GetBool(itf, "-mol2oSubstructure")) mol2oflags|=OEOFlavor::MOL2::Substructure;
  if (GetBool(itf, "-mol2oNameMask"))  mol2oflags|=OEOFlavor::MOL2::NameMask;
  if (GetBool(itf, "-mol2oDefault")) mol2oflags|=OEOFlavor::MOL2::Default;
  if (GetBool(itf, "-smioIsotopes"))   smioflags|=OEOFlavor::SMI::Isotopes;
  if (GetBool(itf, "-smioHydrogens"))  smioflags|=OEOFlavor::SMI::Hydrogens;
  if (GetBool(itf, "-smioRGroups"))    smioflags|=OEOFlavor::SMI::RGroups;
  if (GetBool(itf, "-smioAtomStereo")) smioflags|=OEOFlavor::SMI::AtomStereo;
  if (GetBool(itf, "-smioBondStereo")) smioflags|=OEOFlavor::SMI::BondStereo;
  if (GetBool(itf, "-smioAtomMaps"))   smioflags|=OEOFlavor::SMI::AtomMaps;
  if (GetBool(itf, "-smioCanonical"))  smioflags|=OEOFlavor::SMI::Canonical;
  if (GetBool(itf, "-smioKekule"))     smioflags|=OEOFlavor::SMI::Kekule;
  if (GetBool(itf, "-smioSuperAtoms")) smioflags|=OEOFlavor::SMI::SuperAtoms;
  if (GetBool(itf, "-smioExtBonds"))   smioflags|=OEOFlavor::SMI::ExtBonds;
  if (GetBool(itf, "-smioImpHCount"))  smioflags|=OEOFlavor::SMI::ImpHCount;
  if (GetBool(itf, "-smioSmiMask"))    smioflags|=OEOFlavor::SMI::SmiMask;
  if (GetBool(itf, "-smioDefault")) smioflags|=OEOFlavor::SMI::Default;
  if (GetBool(itf, "-oOEAroModelDaylight"))oflags|=OEOFlavor::Generic::OEAroModelDaylight;
  if (GetBool(itf, "-oOEAroModelOpenEye"))oflags|=OEOFlavor::Generic::OEAroModelOpenEye;
  if (GetBool(itf, "-oOEAroModelTripos")) oflags|=OEOFlavor::Generic::OEAroModelTripos;
  if (GetBool(itf, "-oOEAroModelMMFF")) oflags|=OEOFlavor::Generic::OEAroModelMMFF;
  if (GetBool(itf, "-oOEAroModelMDL")) oflags|=OEOFlavor::Generic::OEAroModelMDL;
  if (GetBool(itf, "-oAroMask")) oflags|=OEOFlavor::Generic::AroMask;
  if (GetBool(itf, "-oRings")) oflags|=OEOFlavor::Generic::Rings;
  if (GetBool(itf, "-oGenericMask")) oflags|=OEOFlavor::Generic::GenericMask;

  ////////////////////////////////////////////
  /// output-flavoring
  ////////////////////////////////////////////
  if (pdboflags)
  {
    if (ofmt!=OEFormat::PDB)
      OEThrow.Fatal("PDB flags incompatible w/ %s",OEGetFormatString(ofmt));
    oflags|=pdboflags;
  }
  if (mdloflags)
  {
    if (ofmt!=OEFormat::MDL && ofmt!=OEFormat::SDF)
      OEThrow.Fatal("MDL flags incompatible w/ %s",OEGetFormatString(ofmt));
    oflags|=mdloflags;
  }
  if (smioflags)
  {
    if (ofmt!=OEFormat::SMI && ofmt!=OEFormat::CAN && ofmt!=OEFormat::ISM && ofmt!=OEFormat::USM)
      OEThrow.Fatal("SMI flags incompatible w/ %s",OEGetFormatString(ofmt));
    oflags|=smioflags;
  }
  if (mol2oflags)
  {
    if (ofmt!=OEFormat::MOL2 && ofmt!=OEFormat::MOL2H)
      OEThrow.Fatal("MOL2 flags incompatible w/ %s",OEGetFormatString(ofmt));
    oflags|=mol2oflags;
  }
  if (mfoflags)
  {
    if (ofmt!=OEFormat::MF)
      OEThrow.Fatal("MF flags incompatible w/ %s",OEGetFormatString(ofmt));
    oflags|=mfoflags;
  }
  if (mmodoflags)
  {
    if (ofmt!=OEFormat::MMOD)
      OEThrow.Fatal("MMOD flags incompatible w/ %s",OEGetFormatString(ofmt));
    oflags|=mmodoflags;
  }
  if (mopacoflags)
  {
    if (ofmt!=OEFormat::MOPAC)
      OEThrow.Fatal("MOPAC flags incompatible w/ %s",OEGetFormatString(ofmt));
    oflags|=mopacoflags;
  }

  if (GetBool(itf, "-oFlavorNone"))
  {
    if (oflags)
      OEThrow.Fatal("-oFlavorNone incompatible with any other output flavors.");

    ofs.SetFlavor(ofmt, 0);
  }
  else if (oflags)
    ofs.SetFlavor(ofmt, oflags|OEOFlavor::Generic::Default);

  if (!OEIsWriteable(ofmt))
    OEThrow.Fatal("non-writeable output format: %s (use -helpformats for supported list)",
                  OEGetFormatString(ofmt));
  else if (!OEIsWriteable(ofmt,ofs.GetFlavor(ofmt)))
    OEThrow.Fatal("writeable output format (%s) but non-writeable flavor: %04X (appropriate -{fmt}oDefault flag may be needed)",
                  OEGetFormatString(ofmt),oflags);

  return true;
}

///////////////////
// Helper functions
string GetInputFile(OEInterface &itf)
{
  if (itf.Has<string>("-ifmt"))
  {
    unsigned int ifmt=OEGetFileType(itf.Get<string>("-ifmt").c_str());
    if (ifmt==OEFormat::UNDEFINED)
      OEThrow.Fatal("Cannot parse -ifmt: %s",itf.Get<string>("-ifmt").c_str());

    string ofname = itf.Get<string>("-ifmt");
    if (GetBool(itf, "-igz"))
      ofname += ".gz";

    return ofname;
  }
  else if (itf.Has<string>("-in"))
    return itf.Get<string>("-in");

  OEThrow.Fatal("No input file specified!");
  return "";
}

unsigned int GetInputFormat(OEInterface &itf)
{
  return OEGetFileType(OEGetFileExtension(GetInputFile(itf).c_str()));
}

//////////////////////////////////////////////////////////
// These are used inside the OpenFormattedIStream template
void OEBKeepROC(oemolstreambase &ifs, OEInterface &itf, unsigned int ifmt)
{
  if (ifmt != OEFormat::OEB || GetOutputFormat(itf) != OEFormat::OEB)
    OEThrow.Fatal("-oeb_keep_roc requires OEB input and output.");

  OEBinaryIOHandlerBase &ihand = ifs.GetBinaryIOHandler();
  ihand.Clear();
  OEInitHandler(ihand,OEBRotCompressOpts(),OEBRotCompressOpts());
}

void MCMolInputFlags(oemolistream &ifs, OEInterface &itf)
{
  if (GetBool(itf, "-mc"))
    ifs.SetConfTest(OEAbsoluteConfTest(GetBool(itf, "-mc_titles")));
  else if (GetBool(itf, "-mc_isomer"))
    ifs.SetConfTest(OEIsomericConfTest(GetBool(itf, "-mc_titles")));

  unsigned int ifmt = GetInputFormat(itf);
  if ((ifmt==OEFormat::XYZ ||
       ifmt==OEFormat::PDB) &&
      (GetBool(itf, "-mc") ||
       GetBool(itf, "-mc_isomer")))
    OEThrow.Warning("PDB, XYZ not multiconformer-capable.");

  if (GetBool(itf, "-oeb_keep_roc"))
    OEBKeepROC(ifs, itf, ifmt);
}

void MCMolInputFlags(oemolithread &ifs, OEInterface &itf)
{
  oeAssert(!GetBool(itf, "-mc"));
  oeAssert(!GetBool(itf, "-mc_isomer"));

  unsigned int ifmt = GetInputFormat(itf);
  if (GetBool(itf, "-oeb_keep_roc"))
    OEBKeepROC(ifs, itf, ifmt);
}

//////////////////////////////////////////////////////////
/// OpenFormattedIStream - open and format the input stream
//////////////////////////////////////////////////////////
template<class StreamType>
void OpenFormattedIStream(StreamType &ifs, OEInterface &itf,
                          string &ifname, OEBabel *babel=0)
{
  unsigned int mol2iflags=0;
  unsigned int pdbiflags=0;
  unsigned int smiiflags=0;
  unsigned int xyziflags=0;
  unsigned int mmodiflags=0;
  unsigned int iflags=0;

  if (!ifs.open(ifname))
    OEThrow.Fatal("Unable to open %s for reading", ifname.c_str());

  unsigned int ifmt = ifs.GetFormat();

  if (babel && ifmt == OEFormat::OEB)
    OEReadHeader(ifs, babel->GetHeader());

  if (GetBool(itf, "-pdbiTER"))       pdbiflags|=OEIFlavor::PDB::TER;
  if (GetBool(itf, "-pdbiEND"))       pdbiflags|=OEIFlavor::PDB::END;
  if (GetBool(itf, "-pdbiENDM"))      pdbiflags|=OEIFlavor::PDB::ENDM;
  if (GetBool(itf, "-pdbiTerMask"))   pdbiflags|=OEIFlavor::PDB::TerMask;
  if (GetBool(itf, "-pdbiALL"))       pdbiflags|=OEIFlavor::PDB::ALL;
  if (GetBool(itf, "-pdbiDATA"))      pdbiflags|=OEIFlavor::PDB::DATA;
  if (GetBool(itf, "-pdbiCHARGE"))    pdbiflags|=OEIFlavor::PDB::CHARGE;
  if (GetBool(itf, "-pdbiRADIUS"))    pdbiflags|=OEIFlavor::PDB::RADIUS;
  if (GetBool(itf, "-pdbiDELPHI"))    pdbiflags|=OEIFlavor::PDB::DELPHI;
  if (GetBool(itf, "-pdbiBasicMask")) pdbiflags|=OEIFlavor::PDB::BasicMask;
  if (GetBool(itf, "-pdbiFormalCrg")) pdbiflags|=OEIFlavor::PDB::FormalCrg;
  if (GetBool(itf, "-pdbiImplicitH")) pdbiflags|=OEIFlavor::PDB::ImplicitH;
  if (GetBool(itf, "-pdbiBondOrder")) pdbiflags|=OEIFlavor::PDB::BondOrder;
  if (GetBool(itf, "-pdbiRings"))     pdbiflags|=OEIFlavor::PDB::Rings;
  if (GetBool(itf, "-pdbiConnect"))   pdbiflags|=OEIFlavor::PDB::Connect;
  if (GetBool(itf, "-pdbiExtraMask")) pdbiflags|=OEIFlavor::PDB::ExtraMask;
  if (GetBool(itf, "-pdbiAllMask"))   pdbiflags|=OEIFlavor::PDB::AllMask;
  if (GetBool(itf, "-pdbiDefault")) pdbiflags|=OEIFlavor::PDB::Default;
  if (GetBool(itf, "-mmodiFormalCrg")) mmodiflags|=OEIFlavor::MMOD::FormalCrg;
  if (GetBool(itf, "-mmodiDefault")) mmodiflags|=OEIFlavor::MMOD::Default;
  if (GetBool(itf, "-mol2iM2H")) mol2iflags|=OEIFlavor::MOL2::M2H;
  if (GetBool(itf, "-mol2iDefault")) mol2iflags|=OEIFlavor::MOL2::Default;
  if (GetBool(itf, "-smiiCanon"))      smiiflags|=OEIFlavor::SMI::Canon;
  if (GetBool(itf, "-smiiStrict"))     smiiflags|=OEIFlavor::SMI::Strict;
  if (GetBool(itf, "-smiiDefault")) smiiflags|=OEIFlavor::SMI::Default;
  if (GetBool(itf, "-xyziFormalCrg")) xyziflags|=OEIFlavor::XYZ::FormalCrg;
  if (GetBool(itf, "-xyziImplicitH")) xyziflags|=OEIFlavor::XYZ::ImplicitH;
  if (GetBool(itf, "-xyziBondOrder")) xyziflags|=OEIFlavor::XYZ::BondOrder;
  if (GetBool(itf, "-xyziRings"))     xyziflags|=OEIFlavor::XYZ::Rings;
  if (GetBool(itf, "-xyziConnect"))   xyziflags|=OEIFlavor::XYZ::Connect;
  if (GetBool(itf, "-xyziExtraMask")) xyziflags|=OEIFlavor::XYZ::ExtraMask;
  if (GetBool(itf, "-xyziDefault")) xyziflags|=OEIFlavor::XYZ::Default;
  if (GetBool(itf, "-iOEAroModelDaylight"))iflags|=OEIFlavor::Generic::OEAroModelDaylight;
  if (GetBool(itf, "-iOEAroModelOpenEye"))iflags|=OEIFlavor::Generic::OEAroModelOpenEye;
  if (GetBool(itf, "-iOEAroModelTripos")) iflags|=OEIFlavor::Generic::OEAroModelTripos;
  if (GetBool(itf, "-iOEAroModelMMFF")) iflags|=OEIFlavor::Generic::OEAroModelMMFF;
  if (GetBool(itf, "-iOEAroModelMDL")) iflags|=OEIFlavor::Generic::OEAroModelMDL;
  if (GetBool(itf, "-iAroMask")) iflags|=OEIFlavor::Generic::AroMask;
  if (GetBool(itf, "-iRings")) iflags|=OEIFlavor::Generic::Rings;
  if (GetBool(itf, "-iGenericMask")) iflags|=OEIFlavor::Generic::GenericMask;

  ////////////////////////////////////////////
  /// input-flavoring
  ////////////////////////////////////////////
  if (mol2iflags)
  {
    if (ifmt!=OEFormat::MOL2)
      OEThrow.Fatal("MOL2 flags incompatible w/ %s",OEGetFormatString(ifmt));
    iflags|=mol2iflags;
  }
  if (pdbiflags)
  {
    if (ifmt!=OEFormat::PDB)
      OEThrow.Fatal("PDB flags incompatible w/ %s",OEGetFormatString(ifmt));
    iflags|=pdbiflags;
  }
  if (smiiflags)
  {
    if (ifmt!=OEFormat::SMI && ifmt!=OEFormat::CAN && ifmt!=OEFormat::ISM &&  ifmt!=OEFormat::USM)
      OEThrow.Fatal("SMI flags incompatible w/ %s",OEGetFormatString(ifmt));
    iflags|=smiiflags;
  }
  if (xyziflags)
  {
    if (ifmt!=OEFormat::XYZ)
      OEThrow.Fatal("XYZ flags incompatible w/ %s",OEGetFormatString(ifmt));
    iflags|=xyziflags;
  }
  if (mmodiflags)
  {
    if (ifmt!=OEFormat::MMOD)
      OEThrow.Fatal("XYZ flags incompatible w/ %s",OEGetFormatString(ifmt));
    iflags|=mmodiflags;
  }

  if (GetBool(itf, "-iFlavorNone"))
  {
    if (iflags)
      OEThrow.Fatal("-iFlavorNone incompatible with any other input flavors.");

    ifs.SetFlavor(ifmt, 0);
  }
  else if (iflags)
    ifs.SetFlavor(ifmt, iflags|OEIFlavor::Generic::Default);

  MCMolInputFlags(ifs, itf);

  if (!OEIsReadable(ifmt))
    OEThrow.Fatal("non-readable input format: %s (use -helpformats for supported list)",
                  OEGetFormatString(ifmt));
  else if (!OEIsReadable(ifmt,ifs.GetFlavor(ifmt)))
    OEThrow.Fatal("readable input format (%s) but non-readable flavor: %04X (appropriate -{fmt}iDefault flag may be needed)",
                  OEGetFormatString(ifmt),iflags);
}

/////////////////////////////////////////////////////////////////////////////
template<class InputStreamType, class OutputStreamType, class MolType, class MolBaseT>
unsigned int ConvertFile(InputStreamType &ifs,
                         OEInterface &itf,
                         OEBabel &babel,
                         string &ofname,
                         unsigned int nmax)
{
  OEDots dots(0,20,"mols"); // Will only count Total by default
  if (OEThrow.GetLevel() <= OEErrorLevel::Verbose)
    dots.SetBigStep(1000);

  bool firstonly=GetBool(itf, "-firstonly");
  unsigned int skip=GetInt(itf, "-skip");
  unsigned int n=GetInt(itf, "-n");
  unsigned int ocount=0;

  OutputStreamType ofs;
  MolType mol;
  bool opened = false;
  if (!nmax && ofname != "")
    opened = OpenFormattedOStream(ofs, itf, ofname, &babel);
  while (OEReadMolecule(ifs,(MolBaseT &)mol))
  {
    if (!opened && ofname != "")
      opened = OpenFormattedOStream(ofs, itf, ofname, &babel);

    dots.Update();

    ICOUNT++;
    if (ICOUNT <= skip)
      continue;

    babel.ApplyTransform((MolBaseT &)mol);

    if (ofname != "")
      ocount += babel.WriteMolecule(ofs, (MolBaseT &)mol);

    if (n > 0 && (ICOUNT-skip) == n)
    {
      OEThrow.Msg(OEErrorLevel::Verbose, "limit reached: %d",n);
      break;
    }
    if (firstonly)
      break;
    if (nmax && nmax <= ocount)
      break;
  }

  if (ocount && OEThrow.GetLevel() <= OEErrorLevel::Verbose)
    dots.Total();

  return ocount;
}

/////////////////////////////////////////////////////////////////////////////
/// Chunker() - split output into several files.
/////////////////////////////////////////////////////////////////////////////
template<class InputStreamType, class OutputStreamType, class MolType, class MolBaseT>
static bool Chunker(OEInterface &itf, OEBabel &babel)
{
  // Determine chunk name prefix and file extension
  string prefix, ext;

  string outpath = GetOutputFile(itf);
  ext = OEGetFileExtension(outpath.c_str());
  prefix = outpath.substr(0,outpath.length()-ext.length()-1);

  if (HasString(itf, "-prefix"))
    prefix = GetString(itf, "-prefix");

  // could be empty because -out is preceded by a .
  if (prefix.empty())
    prefix = "babel_output";

  OEThrow.Msg(OEErrorLevel::Verbose, "chunk prefix: %s", prefix.c_str());
  OEThrow.Msg(OEErrorLevel::Verbose, "chunk file extension: %s", ext.c_str());

  // Determine how many molecules per chunk
  string ifname = GetInputFile(itf);
  unsigned int chunksize=0;
  if (HasInt(itf, "-nchunks"))
  {
    if (HasInt(itf, "-chunksize"))
      OEThrow.Fatal("-chunksize and -nchunks are mutually exclusive");

    unsigned int nchunks = GetInt(itf, "-nchunks");
    if (nchunks < 2)
      OEThrow.Fatal("-nchunks must be more than 1.");

    OEThrow.Info("counting input mols...");
    unsigned int i=0;
    unsigned int n=GetInt(itf, "-n");
    unsigned int skip=GetInt(itf, "-skip");

    if (ifname == OEGetFileExtension(ifname.c_str()))
      OEThrow.Fatal("-nchunks can only be used on actual files");

    MolType mol;
    InputStreamType ifs;
    OpenFormattedIStream(ifs, itf, ifname);
    while (OEReadMolecule(ifs, (MolBaseT &)mol))
      if (++i==skip+n && n)
        break;

    OEThrow.Msg(OEErrorLevel::Debug, "found %u mols", i);
    OEThrow.Msg(OEErrorLevel::Debug, "skipping %u mols", skip);
    i -= skip;

    chunksize = (i/nchunks) + ((i % nchunks)?1:0);
    OEThrow.Info("chunking %d mols", i);
    OEThrow.Info("output chunks will contain %d mols each", chunksize);
  }
  else if (HasInt(itf, "-chunksize"))
    chunksize = GetInt(itf, "-chunksize");

  InputStreamType ifs;
  OpenFormattedIStream(ifs, itf, ifname, &babel);

  unsigned int chunk=0;
  unsigned int ocount=0;
  while (true)
  {
    // Construct chunk filename
    string ofname = prefix + '_' + OENumberToString(chunk+1) + '.' + ext;

    // Do the conversion
    unsigned int j=ConvertFile<InputStreamType,OutputStreamType,MolType,MolBaseT>
                              (ifs,itf,babel,ofname,chunksize);
    ocount += j;

    // Report some stats
    if (j)
    {
      chunk++;
      OEThrow.Info("mols written to %s: %d",ofname.c_str(),j);
    }

    // Determine whether we're done
    if (j < chunksize)
      break;

    unsigned int n = GetInt(itf, "-n");
    if (n && n <= ocount)
      break;
  }
  OEThrow.Info("files created: %d",chunk);
  return true;
}

/////////////////////////////////////////////////////////////////////////////
static void PrintBanner()
{

  OEThrow.Info("          :jGf:           ");
  OEThrow.Info("        :jGDDDDf:         ");
  OEThrow.Info("      ,fDDDGjLDDDf,         ______       _          _");
  OEThrow.Info("    ,fDDLt:   :iLDDL;       | ___ \\     | |        | |");
  OEThrow.Info("  ;fDLt:         :tfDG;     | |_/ / __ _| |__   ___| |");
  OEThrow.Info(",jft:   ,ijfffji,   :iff    | ___ \\/ _` | '_ \\ / _ \\ |");
  OEThrow.Info("     .jGDDDDDDDDDGt.        | |_/ / (_| | |_) |  __/ |");
  OEThrow.Info("    ;GDDGt:''':tDDDG,       \\____/ \\__,_|_.__/ \\___|_|");
  OEThrow.Info("   .DDDG:       :GDDG.    ");
  OEThrow.Info("   ;DDDj         tDDDi      %s - molecular structure file conversion",PROGNAME.c_str());
  OEThrow.Info("   ,DDDf         fDDD,      version: %s",BABEL_VERSION.c_str());
  OEThrow.Info("    LDDDt.     .fDDDj     ");
  OEThrow.Info("    .tDDDDfjtjfDDDGt        Copyright (C) 2005,2006,2007,2008");
  OEThrow.Info("      :ifGDDDDDGfi.         OpenEye Scientific Software, Inc.");
  OEThrow.Info("          .:::.             OEChem version: %s",OEChemGetRelease());
  OEThrow.Info("  ......................    platform: %s",OEChemGetPlatform());
  OEThrow.Info("  DDDDDDDDDDDDDDDDDDDDDD    built: %d",OEChemGetVersion());
  OEThrow.Info("  DDDDDDDDDDDDDDDDDDDDDD  ");
  OEThrow.Info("");

  string licensee;
  if (OEChemGetLicensee(licensee) && !licensee.empty())
    OEThrow.Info("\tlicensee: %s",licensee.c_str());

  string site;
  if (OEChemGetSite(site) && !site.empty())
    OEThrow.Info("\t    site: %s",site.c_str());

  unsigned int expdate[3];
  if (OEChemIsLicensed(0, expdate))
    if (expdate[2]||expdate[1]||expdate[0])
      OEThrow.Info("\t license: expires %d-%02d-%02d",expdate[2],expdate[1],expdate[0]);

  OEThrow.Info("");
  return;
}

static int PrintSupportedFormats()
{
  OEThrow.Info("\tsupported formats:");
  OEThrow.Info("\text          |code| format                |read? |write?");
  OEThrow.Info("\t-------------+----+-----------------------+------+------");
  for (unsigned int code=1;code<OEFormat::MAXFORMAT;++code)
    OEThrow.Info("\t%-12s | %2d | %-21s | %-4s | %-4s",
                 OEGetFormatExtension(code),
                 code,
                 OEGetFormatString(code),
                 OEIsReadable(code)?"yes":"no",
                 OEIsWriteable(code)?"yes":"no");

  OEThrow.Info("\t-------------+----+-----------------------+------+------");
  return 0;
}

/////////////////////////////////////////////////////////////////////////////
int main(int argc,char *argv[])
{
  OEInterface itf;
  OEConfigure(itf,(unsigned char *)InterfaceData);
  OEPlatform::oeosstream oehelp;
  if (OECheckHelp(itf,argc,argv,true,oehelp))
  {
    PrintBanner();
    OEThrow.Info(oehelp.str());
    return 0;
  }
  if (!OEParseCommandLine(itf,argc,argv))
    OEThrow.Fatal("Unable to parse command line.");

  // Determine the level of verbosity (default to Info and above)
  OEThrow.SetStrict(false);
  OEThrow.SetLevel(OEErrorLevel::Info);
  // Checking the mutual exclusion of these options is more code than it's worth
  if (GetBool(itf, "-nowarn"))
    OEThrow.SetLevel(OEErrorLevel::Error);
  if (GetBool(itf, "-quiet"))
    OEThrow.SetLevel(OEErrorLevel::Warning);
  if (GetBool(itf, "-v"))
    OEThrow.SetLevel(OEErrorLevel::Verbose);
  if (GetBool(itf, "-vv"))
    OEThrow.SetLevel(OEErrorLevel::Debug);

  PrintBanner();

  if (GetBool(itf, "-helpformats"))
    return PrintSupportedFormats();

  if (GetBool(itf, "-add2d"))
    if (!OEDepict::OEDepictIsLicensed())
      OEThrow.Fatal("No license found for \"oedepict\"; required for -add2d.");

  if (!(itf.Has<string>("-in") || itf.Has<string>("-ifmt")))
    OEThrow.Fatal("Input file or format (-in or -ifmt) are required.");

  oeosstream os;
  OEWriteSettings(itf,os,true);
  OEThrow.Msg(OEErrorLevel::Debug, "Using the following settings:");
  OEThrow.Msg(OEErrorLevel::Debug, os.str());

  // Write setup to a file
  if (itf.Has<string>("-output_params"))
  {
    oeofstream ofs(itf.Get<string>("-output_params"));
    ofs << os.str();
  }

  // Get the input and output formats
  unsigned int ifmt = GetInputFormat(itf);
  unsigned int ofmt = GetOutputFormat(itf);

  ///////////////////////////////////////////////
  // Determine functions to apply to the molecule
  ///////////////////////////////////////////////
  OEBabel babel(argc, argv, itf);

  // Hydrogen handling
  string hydrogens=itf.Get<string>("-hydrogens");
  if (hydrogens != "same")
  {
    if (IsMCMolInput(itf, ifmt))
      babel.AddTransform(new OEHydrogens<OEMCMolBase>(hydrogens));
    else
      babel.AddTransform(new OEHydrogens<OEMolBase>(hydrogens));
  }

  // Delete invalid stereo
  if (GetBool(itf, "-delete_invalid_stereo") ||
      GetBool(itf, "-canonical_isomeric_smiles") ||
      GetBool(itf, "-canonical_kekule"))
    babel.AddTransform(new OEDeleteInvalidStereo());

  // Stereo from 3d
  if (GetBool(itf, "-stereofrom3d"))
    babel.AddTransform(new OEStereoFrom3D(ofmt));

  // Correct bond stereo
  if (GetBool(itf, "-mdlcorrectstereo"))
    babel.AddTransform(new OEMDLCorrectStereo(ifmt));

  // Copy sd tag to title
  if (itf.Has<string>("-sd2title"))
  {
    if (IsMCMolInput(itf, ifmt))
      babel.AddTransform(new OESD2TitleMC(itf.Get<string>("-sd2title")));
    else
      babel.AddTransform(new OESD2Title(itf.Get<string>("-sd2title")));
  }

  // Perceive residues
  unsigned int perceive_residues_mask=PerceiveResiduesMask(itf);
  if (GetBool(itf, "-perceive_residues") || perceive_residues_mask)
    babel.AddTransform(new OEPerceiveResiduesFunc(perceive_residues_mask));

  // Canonical kekulize
  if (GetBool(itf, "-canonical_kekule"))
    babel.AddTransform(new OECanonicalKekule());

  // 2D depict the coordinates
  if (GetBool(itf, "-add2d"))
  {
    if (!Is3DFormat(ofmt))
      OEThrow.Fatal("Output format incompatible with -add2d: %s",
                    OEGetFormatString(ofmt));
    babel.AddTransform(new DepictCoordinates(ofmt));
  }

  // Output titles to a file
  if (itf.Has<string>("-output_names"))
    babel.AddTransform(new OEOutputNames(itf.Get<string>("-output_names")));

  // Count conformers based on the input format
  OEConfCountBase *iconfcount=0;
  if (Is3DFormat(ifmt))
  {
    if (IsMCMolInput(itf, ifmt))
    {
      OEConfCountMC *c = new OEConfCountMC();
      babel.AddTransform(c);
      iconfcount = c;
    }
    else
    {
      OEConfCount *c = new OEConfCount();
      babel.AddTransform(c);
      iconfcount = c;
    }
  }


  bool useThreadedIO = true;
  if (!OEHasTokenizer(ifmt) ||
      GetBool(itf, "-singlethread") ||
      GetBool(itf, "-mc") ||
      GetBool(itf, "-mc_isomer"))
    useThreadedIO = false;

  ///////////////////////////////////////
  // Do the conversion inside the Chunker
  ///////////////////////////////////////
  string ifname = GetInputFile(itf);
  string ofname = GetOutputFile(itf);
  if (itf.Has<int>("-nchunks") || itf.Has<int>("-chunksize"))
  {
    if (useThreadedIO)
    {
      if (IsMCMolInput(itf, ifmt))
        Chunker<oemolithread,oemolothread,OEMol,OEMCMolBase>(itf,babel);
      else
        Chunker<oemolithread,oemolothread,OEGraphMol,OEMolBase>(itf,babel);
    }
    else
    {
      if (IsMCMolInput(itf, ifmt))
        Chunker<oemolistream,oemolostream,OEMol,OEMCMolBase>(itf,babel);
      else
        Chunker<oemolistream,oemolostream,OEGraphMol,OEMolBase>(itf,babel);
    }
  }
  else  // Do the actual conversion
  {
    if (useThreadedIO)
    {
      oemolithread ifs;
      OpenFormattedIStream(ifs, itf, ifname, &babel);

      if (IsMCMolInput(itf, ifmt))
        ConvertFile<oemolithread,oemolothread,OEMol,OEMCMolBase>(ifs,itf,babel,ofname,0);
      else
        ConvertFile<oemolithread,oemolothread,OEGraphMol,OEMolBase>(ifs,itf,babel,ofname,0);
    }
    else
    {
      oemolistream ifs;
      OpenFormattedIStream(ifs, itf, ifname, &babel);

      if (IsMCMolInput(itf, ifmt))
        ConvertFile<oemolistream,oemolostream,OEMol,OEMCMolBase>(ifs,itf,babel,ofname,0);
      else
        ConvertFile<oemolistream,oemolostream,OEGraphMol,OEMolBase>(ifs,itf,babel,ofname,0);
    }
  }

  // Print information about the input
  oemolistream ifs;
  string ext = ".";
  ext += OEGetFileExtension(ifname.c_str());
  if (!ifs.open(ext))
    OEThrow.Fatal("Trying to read unsupported format, shouldn't get this far");
  unsigned int defaultFlavor = ifs.GetFlavor(ifs.GetFormat());
  OpenFormattedIStream(ifs, itf, ext);
  unsigned int actualFlavor = ifs.GetFlavor(ifs.GetFormat());

  OEThrow.Info("%s: input format: %s%s%s",
               PROGNAME.c_str(),
               (OEIsGZip(ifname.c_str())?"gzipped ":""),
               OEGetFormatString(ifmt),
               ((actualFlavor != defaultFlavor)?" (flavored)":""));

  OEThrow.Info("%s: input flavor: %04X%s",
               PROGNAME.c_str(),
               actualFlavor,
               ((actualFlavor == defaultFlavor)?" (default)":""));

  OEThrow.Info("%s: mols in: %d",PROGNAME.c_str(),ICOUNT);


  if (GetBool(itf, "-molcount"))
    oeout << ICOUNT << oeendl;

  if (iconfcount)
    OEThrow.Info("%s: confs in: %d",PROGNAME.c_str(),iconfcount->GetCount());

  // Print information about the output
  if (itf.Has<string>("-out") || itf.Has<string>("-ofmt"))
  {
    oemolostream ofs;
    string oext = ".";
    oext += OEGetFileExtension(ofname.c_str());
    if (!ofs.open(oext))
      OEThrow.Fatal("Trying to write unsupported format, shouldn't get this far");

    unsigned int defaultFlavor = ofs.GetFlavor(ofs.GetFormat());
    OpenFormattedOStream(ofs, itf, oext);
    unsigned int actualFlavor = ofs.GetFlavor(ofs.GetFormat());

    OEThrow.Info("%s: output format: %s%s%s",
                 PROGNAME.c_str(),
                 (OEIsGZip(ofname.c_str())?"gzipped ":""),
                 OEGetFormatString(ofmt),
                 ((actualFlavor != defaultFlavor)?" (flavored)":"") );

    OEThrow.Info("%s: output flavor: %04X%s",
                 PROGNAME.c_str(),
                 actualFlavor,
                 ((actualFlavor == defaultFlavor)?" (default)":""));

    OEThrow.Info("%s: mols out: %d",PROGNAME.c_str(),babel.GetOCount());

    if (Is3DFormat(ofmt))
      OEThrow.Info("%s: confs out: %d",PROGNAME.c_str(),babel.GetOConfCount());
  }

  // Whether any problems were encountered
  OEThrow.Info("%s: errors: %d",
               PROGNAME.c_str(),OEThrow.Count(OEErrorLevel::Error));
  OEThrow.Info("%s: warnings: %d",
               PROGNAME.c_str(),OEThrow.Count(OEErrorLevel::Warning));

  return OEThrow.Count(OEErrorLevel::Error);
}
