
void findfilenames(string &gridfname, string &cellfname, string &initvalfname, string &schparamfname, string &logfname, string &datafname, string& statefname){
	ifstream perenos("perenos");
	if(!perenos){
		std::cerr<<"function 'findfilenames' cannot open file perenos\n";
		return;
	}
	int i = 0;
	string s;
	while(getline(perenos, s) != NULL){
		i++;
		if(i == 5){
			s.erase( remove(s.begin(), s.end(), '\r'), s.end() );
			gridfname.assign(s);
		}
		if(i == 6){
			s.erase( remove(s.begin(), s.end(), '\r'), s.end() );
			cellfname.assign(s);
		}
	}
	ifstream grid(gridfname.c_str());
	if(!grid){
		std::cerr<<"function 'findfilenames' grid fails file "<<gridfname<<"\n";
		exit(1);
	}
	grid.close();
	perenos.close();
	initvalfname.assign("initval");
	schparamfname.assign("schparam");
	logfname.assign("log");
	datafname.assign("data");
	statefname.assign("state");
	return;
}
