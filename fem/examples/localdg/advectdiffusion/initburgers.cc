
template <class GridType>
class U0 {
public:
  enum { dimDomain = GridType::dimensionworld };  
  typedef FieldVector<double,dimDomain> DomainType;
  typedef FieldVector<double,1> RangeType;
  U0(double eps,bool diff_timestep=true) :
    velocity(0), epsilon(eps), diff_tstep(diff_timestep) {
		myName = "Burgers-Diffusion";
	}
  void evaluate(const DomainType& arg, RangeType& res) const {
    evaluate(0,arg,res);
  }
  void evaluate(double t,const DomainType& arg, RangeType& res) const {
		if(epsilon < 1e-9)
			res = ((2.0*arg[0]-1.0)<0.0)? 1.0:-1.0;
		else
			res = -tanh((2.0*arg[0]-1.0)/(2.*epsilon));
  }
	
	void printmyInfo(string filename)
	{
  std::ostringstream filestream;
  filestream << filename;

  std::ofstream ofs(filestream.str().c_str(), std::ios::app);
	
	ofs << "Problem: " << myName << "\n\n"
			<< "Exact solution: $u(x,y,z,t):=\\displaystyle{-\\tanh\\left( \\frac{x}{2\\varepsilon} \\right) }$\\\\\n\n";
	
	ofs.close();
	}
	
	
  DomainType velocity;
  double epsilon;
  bool diff_tstep;
	string myName;
};


