#ifndef CPHEXFUNCTIONS5_H
#define CPHEXFUNCTIONS5_H

//#include <eigen3/Eigen/Dense>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <complex>

typedef std::mt19937 MyRNG;             //Mersenne Twister random number gen
uint32_t seed_val;                      //For the generator

unsigned seedv = std::chrono::system_clock::now().time_since_epoch().count();

MyRNG rng;                              // One instance

//We set up the random number generator

void initialize()
{

        rng.seed(seedv);

}

using std::vector;	using std::ifstream;

struct Point {

        int indx;               // The indx will store the grid point
        int cnum;               // cnum will be the cell number
        int ctype;              // ctype will denote the cell type
        int inBulk;             // inBulk 1 if the cell is a bulk grid point
        int Peri;               // Peri 1 if the cell is a perimeter point
	int Vertex;		// Vertex 1 if the cell is a vertex
	int hexrw;		//0 for an even layer and 1 for an odd layer
	double xcart;		//X coordinate in cartesian
	double ycart;		//Y coordinate in cartesian
	double xhex;		//X coordinate in a hex lattice	
	double yhex;		//Y coordinate in a hex lattice
	int layer;		// This keeps track of which layer
	double Activ;		// This keeps track of the activity of the point
	int edgemP;		// This tracks whether the top and bottom layersat the point are both edge
	double zhex;		//We track the distance between layers
	int nperi;		//The number of perimeter interfaces
	int nedge;		// Number of matching edges at a point
	int neiind[6];		// The index of the neighbors
	int neitype[6];		// The type of the corresponding neighbor
	int quadr;		// The quadrant the point lies in
	double Bi;		// Bilayer coupling strength
	double Kapp;		// Curvature coefficient
	double theta2;		// Site diameter will approximate arc length
	void initpts(int lay) 
	{
		cnum = 0; ctype = 0; layer = lay; Activ = 0; edgemP = 0; Bi = 0; 
	}

};

//bookmark
///////////////////////////
//We keep track of grid properties


struct Grid {


        int Lx;                 // The x length of the grid
	int Ly;			// The y length of the grid
        int NP;                 // The total number of grid points	
        int Nnum;               // The total number of cells
        int Ntyp;               // The total number of cell types
	int rows;		// The total number of rows
	int colmn;		// The total number of columns
	int Mag;		// The magnitization for the Ising model
	double Lat;		//The spacing between lattice points
	double Lxhex;              // The x length in hex coordinates
	double Lyhex;		// The y length in hex coordinates 
	double COM;		//The COM of the grid
	double COMx;		//The x coordinate of the COM
	double COMy;		//The y coordinate of the COM
	double Xlo;		//The left boundary
	double Xhi;		//The right boundary
	double Ylo;		//The bottom boundary
	double Yhi;		//The top boundary
	void gdset1(int a, int b, int c) 
	{
		Lx = a; NP = b; Nnum = c; 
		if ( (b %2) == 0 ) { rows = 0;}
	        if ( (b %2) == 1 ) { rows = 1;} 
	} 
	void gdset2(int a, int b, int c, double d) 
	{
		Lx = a; NP = b; Nnum = c; Lat = d;
		if ( (b %2) == 0 ) { rows = 0;}
	        if ( (b %2) == 1 ) { rows = 1;} 
		Lxhex = ((NP % Lx)+1)*Lat*sqrt(3)/2; 
		COM = (2*Lat*(Lx+1)+1)/4;	
		Xlo = 1.5*Lat; Xhi = Lx*Lat;
	}


	void gdset3(int a, int b, int c, double d)
	{
		Lx = a; Ly = b; NP = a*b;  Nnum = c; Lat = d;
		if ( (b %2 ) == 0 ) { rows = 0;}
		if ( (b %2 ) != 0 ) { rows = 1;}
		double rad3 = sqrt(3);
		Xlo = Lat; Xhi = (Lx+0.5)*Lat*rad3/2;
		Ylo = Lat*rad3/2; Yhi = Ly*Lat*rad3/2;
		Lxhex = Lx*Lat; Lyhex = Ly*Lat;
		COM = (Xhi+Xlo)/2;
		COMx = (Xhi+Xlo)/2;
		COMy = (Yhi+Ylo)/2;
	}	

};

/////////////////////////////////////////////////
//We'll also keep track of the points with different neighbors

struct Nlist {

	int n;		//The size of the neighbor list, should be even		
	std::vector<int> List;	//The indices that can change cells		
	void showlist()	//This will show elements of the neighbor list
	{
		int i = 0;
		for (i=0; i<n; i++)
			std::cout << List[i] << std::endl;
	}
};

/////////////////////////////////////////////
//This will keep track of the cells in the bending regimes
//bookmark

struct Bendlist {
	int nreg1;	//The number of cells in region 1
	int nreg2;	//The number of cells in region 2
	int nactive;	//The number of cells with bilayer coupling 
	std::vector<int> reg1list;	//The indices of cells in region 1
	std::vector<int> reg2list;	//The indices of cells in region 2

	void initialize()
	{
		nreg1 = 0; nreg2 = 0; nactive = 0;
		reg1list.clear(); reg2list.clear();
	}
};


//bookmark
/////////////////////////////////////
//We keep track of cell properties
//
struct Cellprop {

        double lam1;            // The scale factor for area
        double Ai;              // The current cell area
	double Ai2;
	double Aihex;		// The area in hex coordinates
	double Aihex2;
        double A0;              // The desired cell area
	double Ainit;		// The initial cell area
        double lam2;            // The scale factor for perimeter
        double Pi;              // The current cell perimeter
	double Pi2;
	double Pihex;		// The perimeter in hex coordinates 
	double Pihex2;
	double Pinit;		// The initial cell perimeter
        double P0;              // The desired cell perimeter
        double type;            // The cell type
        double Ji;              // Adhesion strength
	int lay;		// The grid layer the cell is on
	double Jlay;		// Bilayer adhesion
	double lamact;		//Scaling for activity
	double maxact;		//Max value for activity
	double edgemC;		//The number of matching edges between layers
	double edgemC1;
	double edgemC2;
	double matchratio;	//The ratio of matching edges to edges
	double xcom;		// Center of mass x
	double ycom;		// Center of mass y
	double kapp;		// Scaling factor for curvature
	double PAi;		// Current Shape index
	double PA0;		// Desired Shape index
	int quad;		//quadrant of the com
	int connectx;		//Will keep track of connectivity
	int connecty;
	double disconx;		//Location of x discontinuity
	double discony;		//Location of y discontinuity
	double hexaorder;	// The hexatic order index
	vector<int> neighs;	//The neighboring cell indices

	void getAhex(double L)
	{
                Aihex = Ai*sqrt(3)*L*L/2; 
	}
        void getPhex(double L)
	{
                Pihex = Pi*L/sqrt(3); 
	}
        double Asd()
	{
                return lam1*(Ai-A0)*(Ai-A0); 
	}
	double Asd2()
	{
		return lam1*(Ai2-A0)*(Ai2-A0); 
	}
	double Asdhex()
	{
                return lam1*(Aihex-A0)*(Aihex-A0); 
	}
	double Asdhex2()
	{
                return lam1*(Aihex2-A0)*(Aihex2-A0); 
	}
        double Psd()
	{
                return lam2*(Pi-P0)*(Pi-P0); 
	}
	double Psd2()
	{
		return lam2*(Pi2-P0)*(Pi2-P0); 
	}
	double Psdhex()
	{ 
                return lam2*(Pihex-P0)*(Pihex-P0); 
	}
	double Psdhex2()
	{
                return lam2*(Pihex2-P0)*(Pihex2-P0); 
	}
	void getshapei()
	{
		if (Aihex != 0) { PAi = Pihex/sqrt(Aihex); }
		if (Aihex == 0) { PAi = 0;}  	
	}
	void getshape0()
	{
		if ( A0!= 0 ) { PA0 = P0/sqrt(Aihex);}
		if ( A0 == 0 ) { PA0 = 0; } 	
	}
        void CPclear()
	{
                Ai = 0; Pi = 0; Ai2 = 0; Pi2 = 0; Aihex = 0; Aihex2 = 0; Pihex = 0; Pihex2 = 0; edgemC = 0; edgemC1 = 0; edgemC2 = 0; hexaorder = 0;
	}
	void CPclearscl() 
	{
		lam1 = 0; lam2 = 0; Jlay = 0; 
	}
	void CPinit()
	{
		lam1 = 0; Ai = 0; Aihex = 0; Ai2 = 0; Aihex2 = 0; A0 = 0; lam2 = 0; Pi = 0; Pihex = 0; Pi2 = 0; Pihex2 = 0; P0 = 0; edgemC = 0; kapp = 0; connectx = 0; connecty = 0; disconx = 0; discony = 0; hexaorder = 0; xcom = 0; ycom = 0;
	}
	void CPset(double a, double b, double c, double d)
	{ 
		lam1 = a; A0 = b; lam2 = c; P0 = d; 
	} 

	void CPsetall(double a, double b, double c, double d, double e, double f, double g, double h)
	{ 
		lam1 = a; A0 = b; lam2 = c; P0 = d; Ji = e; Jlay = f; lamact = g; maxact = h; if ( h <= 0 ) { maxact = 1;}
	} 

	void CPsetallK(double a, double b, double c, double d, double e, double f, double g, double h, double j)
        { 
		lam1 = a; A0 = b; lam2 = c; P0 = d; Ji = e; Jlay = f; lamact = g; maxact = h; if (h <= 0) { maxact = 1; } kapp = j;
	}
	void Cellupdate()
	{	
	Ai = Ai2; Pi = Pi2; edgemC = edgemC2; Aihex = Aihex2; Pihex = Pihex2; 
	if (Pi != 0)
		{matchratio = (double)edgemC/(double)Pi; }
	else
		{matchratio = 0;}

	}  
	void Cellupdate2()
	{ 
	Ai2 = Ai; Aihex2 = Aihex; Pi2 = Pi; Pihex2 = Pihex; edgemC2 = edgemC;
	if (Pi2 != 0)
		{matchratio = (double)edgemC2/(double)Pi2; }
        else
		{matchratio = 0;}
	}

	void Getmatchratio()
	{
		if (Pi != 0)
		{matchratio = (double)edgemC/(double)Pi; }
        	else
		{matchratio = 0;}

        }


	void Addneigh(int ne)
	{ 
		
		int chk = 0;
		
		if (neighs.empty() )
		{
			neighs.push_back(ne);
		}

		else
		{
			int i = 0, nsze = neighs.size();
			for (i=0; i < nsze; i++ )
			{
				if (neighs[i] == ne)
				{
					chk = 1;
				}
			}

	
		}

		if (chk == 0)
		{
			neighs.push_back(ne);
		}

	}		//End add neighbor


};


//bookmark
///////////////////////////////
//We keep track of the system energy
//
// And this will keep track of the system energies
//
struct Hamil {

        double Hadhold;
        double Hadhnew;

	double Hlayold;
	double Hlaynew;

        double HAold;
        double HAnew;

        double HPold;
	double HPnew;

	double Hactold;
	double Hactnew;

        double Hold;
        double Hnew;
	double dH;

        double getHold() {
		Hold = Hadhold + HAold + HPold + Hlayold + Hactold;
                return Hadhold + HAold + HPold + Hlayold + Hactold; }

	double getHold2() {
		Hold = Hadhold + Hlayold + HAold + HPold + Hactold;
		return Hadhold + Hlayold + HAold + HPold + Hactold; }


        double getHnew() {
		Hnew = Hadhnew + HAnew + HPnew + Hlaynew + Hactnew;
                return Hadhnew + HAnew + HPnew + Hlaynew + Hactnew; }

	double getHnew2() {
		Hnew = Hadhnew + Hlaynew + HAnew + HPnew + Hactnew;
                return Hadhnew + Hlaynew + HAnew + HPnew + Hactnew; }


	void Hupdate() {
		Hold = Hnew; }

	void Hclear() {
        Hadhold = 0; Hlayold = 0; HAold = 0; HPold = 0; Hactold = 0; Hold = 0;
        Hadhnew = 0; Hlaynew = 0; HAnew = 0; HPnew = 0; Hactnew = 0; Hnew = 0; 
	dH = 0;}


};


//////////////////////////////////////////////////////////////////////////////
//Here we get the parameters from file
void Getfileparams(std::string& file, int& npx, int& lat, int& numcells, int& confluent, double epsplane[], double epsbi[], double Aconstrt[], double Pconstrt[], double Actconstrt[], double Actlen[], double& Temp, int& steps, double& Areal1, double& Areal2, double& Peril1, double& Peril2, double& Bilast, double& Bfreq )
{
	std::string words;
	int skip;
	std::ifstream fi(file.c_str());

	fi >> words >> lat >> words >> confluent >> words >> numcells >> words >> skip  >> words >> skip;

	fi >> words >> epsplane[0] >> epsplane[1];
	fi >> words >> epsbi[0] >> epsbi[1];
	fi >> words >> Aconstrt[0] >> Aconstrt[1];
	fi >> words >> Pconstrt[0] >> Pconstrt[1];
	fi >> words >> Actconstrt[0] >> Actconstrt[1];
	fi >> words >> Actlen[0] >> Actlen[1];
	fi >> words >> Temp;
	fi >> words >> steps;
	fi >> words >> Areal1 >> words >> Areal2;
	fi >> words >> Peril1 >> words >> Peril2;





}

///////////////////////////////////////////////////////////////////////////
//Here we get the data from the region file
void Getregsystemsetup(std::string& file, int& npoints, int& npx, int& latck, int& numcells, int& confluent, Point po1[], Point po2[] )
{

	std::string words;
	std::ifstream fi(file.c_str());
	fi >> words >> words >> words >> words >> words;
	fi >> npoints >> npx >> latck >> numcells >> confluent;
	fi >> words >> words >> words >> words >> words;

	int ploop = 0, index, xcoord, ycoord, Celnum, hexrow, skip;
	double Xhex, Yhex, rad3 = std::sqrt(3);

	po1 = new Point[npoints]; po2 = new Point[npoints];

	for (ploop = 0; ploop < npoints; ploop++)
	{
		xcoord = (ploop%npx)+1; ycoord = std::floor(ploop/npx)+1;
		fi >> index; po1[ploop].indx = index; po2[ploop].indx = index;
		po1[ploop].xcart = xcoord; po2[ploop].xcart = xcoord;
		po2[ploop].ycart = ycoord; po2[ploop].ycart = ycoord;
		fi >> Xhex; po1[ploop].xhex = Xhex; po2[ploop].yhex = Yhex;
		fi >> Yhex; po1[ploop].yhex = Yhex; po2[ploop].yhex = Yhex;
		hexrow = ycoord%2;

		if (hexrow==1) {po1[ploop].hexrw = 1; po2[ploop].hexrw = 1; }
		if (hexrow==0) {po1[ploop].hexrw = 1; po2[ploop].hexrw = 1; }

		fi >> po1[ploop].cnum; fi >> po2[ploop].cnum;
		po1[ploop].zhex = 0;	po2[ploop].zhex = 0;

	}

		



}



///////////////////////////////////////////////////////////////////////////
//Here we get data from a dat file
void Getdatsystemsetup(std::string& file,int& npoints,Point po1[],Point po2[], double Shift)
{
	std::string words;
        std::ifstream fi(file.c_str());

	int inum = 0, ids, cnum, peri, bimatch;
	double xh, yh, zh, vx, vy, vz, activ;
        fi >> words >> words >> words >> words >> words >> words >> words >> words >> words >> words >> words;

	for (inum=0; inum < 2*npoints; inum++)
	{

		fi >> ids >> cnum >> xh >> yh >> zh >> vx >> vy >> vz >> peri >> bimatch >> activ; 

		if (inum < npoints)
		{
			po1[inum].indx = ids-1; po1[inum].cnum = cnum;
			po1[inum].xhex = xh; po1[inum].yhex = yh;
			po1[inum].zhex = zh; po1[inum].nperi = peri;
			po1[inum].nedge = bimatch; po1[inum].Activ = activ;

		}

		if (inum >= npoints)
                {
                        po2[inum-npoints].indx = ids-1; po2[inum-npoints].cnum = cnum;
                        po2[inum-npoints].xhex = xh - Shift; po2[inum-npoints].yhex = yh;
                        po2[inum-npoints].zhex = zh; po2[inum-npoints].nperi = peri;
                        po2[inum-npoints].nedge = bimatch; po2[inum-npoints].Activ = activ;

                }
	

	}

	fi.clear();
	fi.close();

}


//////////////////////////////////////
//This will copy the properties of one point and another
void copypt(Point Poit[], Point Poit2[], int a, int b)

{

	Poit[a].cnum = Poit2[b].cnum;

	Poit[a].Peri = Poit2[b].Peri;

 	Poit[a].edgemP = Poit2[b].edgemP;

	Poit[a].zhex = Poit2[b].zhex;

	Poit[a].nedge = Poit2[b].nedge;


}


//////////////////////////////////////
//This will get the row number of each point

void gethxrwall(Point Poit[], Grid gee)
{


int xx = gee.Lx;

int ss = gee.NP;

int hxchk;

int i = 0;
int y;


for (i=0; i<ss; i++)
{

	y = floor(i/xx) + 1;
	hxchk = y %2;

	if (hxchk == 1) { Poit[i].hexrw = 1;} //Odd rows are 1

	if (hxchk == 0) { Poit[i].hexrw = 0; } // Even rows are 0


}


}

////////////////////////////////
///////////////////////////////////////
//bookmark
//We set the neighbor indices points using a hex lattice

vector<int> neighlistset(Point poit[], Grid gee, Cellprop CC[], int a)
{

int tyl = gee.Nnum, xx = gee.Lx, nce[tyl] = {0}, ss = gee.NP, remchk;

int rn,ct = poit[a].cnum, pnum = 0, index, y, hxchk;

remchk = (a+1) % xx;

y = floor(a/xx) + 1;
hxchk = (y%2);

if (hxchk == 1) { poit[a].hexrw = 1;} //Odd rows are 1

if (hxchk == 0) { poit[a].hexrw = 0; } // Even rows are 0

int row = poit[a].hexrw;

//poit[a].Neigh.clear();

vector<int> neighind;
neighind.clear();

//std::cout<< ss << " " << xx << " " << a << ", row = " << row << " " << remchk << std::endl;

//The left
//


if (remchk != 1 ) { index = a-1;}

if (remchk == 1 ) {index = a+xx-1;}


rn = poit[index].cnum; nce[rn] = nce[rn] + 1; neighind.push_back(index);
if (rn != ct) {pnum  = pnum+1; CC[ct].Addneigh(rn);}
poit[a].neiind[0] = index;
poit[a].neitype[0] = rn;

//std::cout<< "left " << index << std::endl;

///////////////////////////////////

//The right
//

if (remchk != 0 ) {index = a+1; }

if (remchk == 0 ) {index = a-xx+1; }

rn = poit[index].cnum; nce[rn] = nce[rn] + 1; neighind.push_back(index);
if (rn != ct) {pnum = pnum+1; CC[ct].Addneigh(rn);}
poit[a].neiind[1] = index;
poit[a].neitype[1] = rn;


//std::cout<<"right " << index << std::endl;

////////////////////////////////////
//Now we start the hard cases
//For a hex lattice, we need to keep track of the row 

////////////////

//Upper Left
//


if  (row == 1 ) //Odd row

{


        if ( (a > xx-1) && (remchk != 1 ) ) {index = a-xx-1; }

// Upper Left of Top Row, 

        if ( (a < xx ) && ( a != 0 ) ) {index = ss-xx+a-1; }

// Upper Left of Left Column
 //
        if ( (remchk == 1 ) && (a != 0 ) ) {index = a-1; }

// Upper left of 0
//
        if ( a == 0 ) { index = ss-1; }

rn = poit[index].cnum; nce[rn] = nce[rn] + 1; neighind.push_back(index);
if (rn != ct) {pnum = pnum + 1;	CC[ct].Addneigh(rn);} 

poit[a].neiind[2] = index;
poit[a].neitype[2] = rn;


}	//End odd row

///////////

if  (row == 0 )		//Even row, the upper left is the index above

{

	if (a > xx-1 ) {index = a-xx; }


	if (a < xx ) {index = ss-xx-a; }

rn = poit[index].cnum; nce[rn] = nce[rn] + 1; neighind.push_back(index);
if (rn != ct) {pnum = pnum+1; CC[ct].Addneigh(rn);}
poit[a].neiind[2] = index;
poit[a].neitype[2] = rn;


}


//std::cout<<"upper left " << index << std::endl;

/////////////////////////////////////////////////
//Upper right


// Upper Right of odd row, is the index above
//

if (row == 1 )

{


	if (a > xx-1 ) {index = a-xx; }


	if (a < xx ) {index = ss-xx+a; }

rn = poit[index].cnum; nce[rn] = nce[rn] + 1; neighind.push_back(index);
if (rn != ct) {pnum = pnum+1; CC[ct].Addneigh(rn);}
poit[a].neiind[3] = index;
poit[a].neitype[3] = rn;



}	//End odd case

////////////

if (row == 0) //Upper right of even row is upper right

{

	if ( (a > xx-1) && (remchk != 0 ) ) { index = a-xx+1; }

// Upper Right of top row
//
        if ( (a < xx) && ( a != xx-1) ) {index = ss-xx+a+1; }
      
 // Upper Right of Right Column
 //
        if ( (remchk == 0 ) && (a != xx-1 ) ) { index = a-2*xx+1; }
 

// Upper Right of xx-1 
//
        if ( a == xx-1 ) { index = ss; }

rn = poit[index].cnum; nce[rn] = nce[rn] + 1; neighind.push_back(index);
if (rn != ct) {pnum = pnum + 1; CC[ct].Addneigh(rn);}

poit[a].neiind[3] = index;
poit[a].neitype[3] = rn;


}

//std::cout<<"upper right " << index << std::endl;

///////////////////////////////////////////////

// Lower Left
//

//
if (row == 1) //Lower left of odd row

{

        if ( (a < ss-xx ) && ( remchk != 1 ) ) { index = a+xx-1; }


// Lower Left of Bottom Edge, need asymmetry for a hex lattice

//        if ( (a > ss-xx ) && ( a != ss-xx ) ) { index = a+xx-ss-1; }

	if ( (a > ss-xx ) && ( a != ss-xx ) ) { index = a+xx-ss; }

// Lower Left of Left Edge
//
        if ( (a < ss-xx ) && (remchk == 1 ) ) { index = a+2*xx-1; }


// Lower Left of left corner
//
//        if ( a == ss-xx ) {index = xx-1; }

	if (a == ss-xx ) { index = 0; }

rn = poit[index].cnum; nce[rn] = nce[rn] + 1; neighind.push_back(index);
if (rn != ct) {pnum = pnum + 1; CC[ct].Addneigh(rn);} 

poit[a].neiind[4] = index;
poit[a].neitype[4] = rn;


}	//End odd case

///////////


if (row == 0) //Lower left of even row, becomes the index below

{
	if (a < ss-xx) { index = a+xx; }

	if (a > ss-xx ) { index = a-ss+xx; }

	if (a == ss-xx ) { index = 0; }

rn = poit[index].cnum; nce[rn] = nce[rn] + 1; neighind.push_back(index);
if (rn != ct) {pnum = pnum+1; CC[ct].Addneigh(rn);}

poit[a].neiind[4] = index;
poit[a].neitype[4] = rn;


}

//std::cout<<"lower left " << index << std::endl;

//////////////////////////////////////////////////


// Lower Right
//
//

if (row == 1 )	//Odd row, lower right becomes index below
{

	// Bottom edge needs special treatment
	
	//if ((a > ss-xx ) && ( a != ss-1)  ) {index = a-ss+xx; }
	if ((a > ss-xx-1) && (a != ss-1 ) ) { index = a-ss+xx+1; }  

	if ( a == ss-1 ) { index = 0; }

	if (a < ss-xx ) { index = a+xx; }

rn = poit[index].cnum; nce[rn] = nce[rn] + 1; neighind.push_back(index);
if (rn != ct) {pnum = pnum+1; CC[ct].Addneigh(rn);}

poit[a].neiind[5] = index;
poit[a].neitype[5] = rn;

}	//End odd case

/////////
//
//

if (row == 0) //Even row

{

        if ( (a < ss-xx ) && ( remchk != 0 ) ) {index = a+xx+1; }

// Lower Right of Bottom Edge

        if ( (a > ss-xx-1 ) && ( a != ss-1 ) ) { index = a-ss+xx+1; }

// Lower Right of Right Edge
//
        if ( (a < ss-xx) && (remchk == 0 ) ) { index = a+1; } 

	// Lower Right of right corner
//
        if ( a == ss-1 ) {index = 0; }


rn = poit[index].cnum; nce[rn] = nce[rn] + 1; neighind.push_back(index);
if (rn != ct) {pnum = pnum + 1; CC[ct].Addneigh(rn);} 

poit[a].neiind[5] = index;
poit[a].neitype[5] = rn;


} //End even row

//std::cout<<"lower right " << index << std::endl;


/////////////////////////////////////


	//We'll add the number of different cell interfaces to the end of the vector containing the neighbor indices for the sake of convienence
	//
//	neighind.push_back(pnum+pnum2);
	neighind.push_back(pnum);



if (pnum > 1) { poit[a].Vertex = 1; poit[a].Peri = 1;  }
if (pnum == 1 ) { poit[a].Peri = 1; }
else {poit[a].Peri = 0; }

poit[a].nperi = pnum;

	return neighind;

}

///////////////////////////////////////////////////
//Once the neighbor list is set we can check the arrays

vector<int> getneighs(Point poit[], int a)
{

	vector<int> neighind;
	neighind.clear();

	int j = 0;

	int c1 = poit[a].cnum;
	int c2;
	int mat = 0;

	for (j=0; j<6; j++)
	{
		c2 = poit[a].neiind[j];
		neighind.push_back(c2);

		if (c1 != c2) {mat = mat+1; }

	}


	neighind.push_back(mat);

	return neighind;

}




//////////////////////////////////////////////////
//This will be a quick way to update site properties

void neighpropset(Point poit[], Point poit2[], Cellprop CC[], Cellprop LL[], int a)
{

	int j = 0;

	int index;

	int cn1 = poit[a].cnum;

	int cn2 = poit2[a].cnum;
       
	int cn3, cn4;

	

//	int per = 0;

	poit[a].nedge = 0;
	poit2[a].nedge = 0;

	int nc = 0;

	for (j = 0; j<6; j++)
	{

	 index = poit[a].neiind[j];

	cn3 = poit[index].cnum;

	cn4 = poit2[index].cnum;

		if (cn1 != cn3 )

		{

			CC[cn1].Addneigh(cn3);

			nc +=1;

			if (cn2 != cn4) 
				{
					poit[a].nedge += 1;
					poit2[a].nedge += 1;
					LL[cn2].Addneigh(cn4);
				}		
		
		}

	}


	if (nc > 0 ) { poit[a].Peri = 1; }

	if (nc == 0 ) { poit[a].Peri = 0; }

	poit[a].nperi = nc;

	poit[a].zhex = poit[a].nedge*CC[cn1].kapp;

	poit2[a].zhex = poit2[a].nedge*LL[cn2].kapp;


}

//////////////////////////////////////////////////
//This will get site properties within a specific region

void neighpropget(Point poit[], Point poit2[], Cellprop CC[], Cellprop LL[], int a, int& nper, int& nmatc)
{

	int j = 0;

	int index;

	int cn1 = poit[a].cnum;

	int cn2 = poit2[a].cnum;
       
	int cn3, cn4;

	

//	int per = 0;

//	poit[a].nedge = 0;
//	poit2[a].nedge = 0;

	int nc = 0;
	int ncm = 0;

	for (j = 0; j<6; j++)
	{

		index = poit[a].neiind[j];

		cn3 = poit[index].cnum;

		cn4 = poit2[index].cnum;


		if (cn1 != cn3 )

		{

//			CC[cn1].Addneigh(cn3);

			nc +=1;

			if (cn2 != cn4) 
				{
					//poit[a].nedge += 1;
					//poit2[a].nedge += 1;
					ncm +=1;
//					LL[cn2].Addneigh(cn4);
				}		
		
		}

	}


//	if (nc > 0 ) { poit[a].Peri = 1; }

//	if (nc == 0 ) { poit[a].Peri = 0; }

	nper += nc;
	nmatc += ncm;

//	poit[a].nperi += nc;

//	poit[a].zhex = poit[a].nedge*CC[cn1].kapp;

//	poit2[a].zhex = poit2[a].nedge*LL[cn2].kapp;


}




//////////////////////////////////////
//We get the coordinates of the point
//

void getcoord(Point pt[], Grid ge, int a)
{

//int tyl = ge.Nnum;
int xx = ge.Lx;
//int ss = ge.NP;

double x, y;

x = (a % xx) + 1;
y = floor(a/xx) + 1;

pt[a].xcart = x;
pt[a].ycart = y;

}


//////////////////////////////

//We convert Cartesion coordinates to those in a hex lattice
//If the distance between hex centers is L, the side length is L/sqrt(3)
//If we set the side length to 1, the distance between centers is sqrt(3)
void getcoordhex(Point pt[], Grid Ge, int a, double L)
{

	double xhexx, yhexx;
//	int tyl = Ge.Nnum;
	int xx = Ge.Lx;
//	int ss = Ge.NP;
	double mdptx = Ge.COMx;
	double mdpty = Ge.COMy;
	int x = (a % xx) + 1;
	int y = floor(a/xx) + 1;
	int hxchk = y %2;

	if (hxchk == 1) { pt[a].hexrw = 1;} //Odd rows are 1

	if (hxchk == 0) { pt[a].hexrw = 0; } // Even rows are 0


	int rown = pt[a].hexrw;

	if (rown == 1)	//Odd rows
	{
		xhexx = (double)x*L;
		yhexx = (double)y*L*sqrt(3)/2;

	}

	if (rown == 0) // Even rows
	{
		xhexx = ((double)x+0.5)*L;
		yhexx = (double)y*L*sqrt(3)/2;
	}

	pt[a].xhex = xhexx; pt[a].yhex =  yhexx;

	if ( (xhexx > mdptx) && ( yhexx > mdpty) ) { pt[a].quadr = 1; }

	if ( (xhexx < mdptx) && ( yhexx > mdpty) ) { pt[a].quadr = 2; }

	if ( (xhexx < mdptx) && ( yhexx < mdpty) ) { pt[a].quadr = 3; }

	if ( (xhexx > mdptx) && ( yhexx < mdpty) ) { pt[a].quadr = 4; }


}

///////////////////
vector<double> index2hex(Grid ge, int a, double L)
{


//	int tyl = ge.Nnum;
	int xx = ge.Lx;
//	int ss = ge.NP;

	int x, y;

	x = (a % xx) + 1;
	y = floor(a/xx) + 1;

	double xhexx, yhexx;
	int rown = y%2;

	vector<double> hxcor;

	hxcor.clear();

	if (rown == 1)	//Odd rows
	{
		xhexx = (double)x*L;
		yhexx = (double)y*L*sqrt(3)/2;

	}

	if (rown == 0) // Even rows
	{
		xhexx = ((double)x+0.5)*L;
		yhexx = (double)y*L*sqrt(3)/2;
	}


	hxcor.push_back(xhexx); hxcor.push_back(yhexx);

	return hxcor;

}


/////////////////////////



//Here we get area using a hex lattice with spacing L

void getAreahex( Cellprop CC[], int k, double L )
{

	CC[k].Aihex = CC[k].Ai*sqrt(3)*L*L/2;

}


//////////////////////////
//

void getPermhex( Cellprop CC[], int k, double L )
{

        CC[k].Pihex = CC[k].Pi*L/sqrt(3);

}

/////////////////////////////////////////////////


///////////////////////////////////////////////
//Gets the bilayer energy on a hex lattice


double layerpermhex(Point poit[], Cellprop CCC[], Point poit2[], Cellprop PPP[], int a, vector<int> n1, vector<int> n2)

{

int i, nind, nind2,  cn1, cn2, cn3, cn4;

double Hlay = 0,J1,J12lay;

vector<int> layneigh;

//We only look at grid points who have neighbors in different cells
//
//nmatch = poit[a].Neigh[8];

cn1 = poit[a].cnum;	//Cell number of first grid point

cn3 = poit2[a].cnum;             //Cell number of second grid point

J1 = CCC[cn1].Jlay;

J12lay = (J1 + PPP[cn3].Jlay)/2;


//We loop over the left, right, top, and bottom neighbors
for (i=0; i<6; i++)

	{

	nind = n1[i];

	cn2 = poit[nind].cnum;

        nind2 = n2[i];        //We get the index of the neighbors
        cn4 = poit2[nind2].cnum;          //We get the cell to which it belongs

//We treat them as edges if the the neighbors are different
//And add the info on a site and cell level
//
	if ( (cn1 != cn2) && (cn3 != cn4) )
	{ Hlay = Hlay + J12lay;}


	} // End the loop over neigbors



//	Hlay = Hlay + J1*nmatch;

	//Now we get the energy between layers

	return Hlay;

}	//End the function

/////////////////////////////////////////////////////

///////////////////////////////
//
////This will track the number of edges

int bilayeredgecnthex(Point poit[], Point poit2[], int a, vector<int> n1, vector<int> n2)

{


int i, nind, nind2, getmatch=0, cn1, cn2, cn3, cn4;

vector<int> layneigh;

//We only look at grid points who have neighbors in different cells
//

cn1 = poit[a].cnum;	//Cell number of first grid point

cn3 = poit2[a].cnum;             //Cell number of second grid point

//We loop over the hex neighbors
for (i=0; i<6; i++)

	{

	nind = n1[i];

	cn2 = poit[nind].cnum;

        nind2 = n2[i];        //We get the index of the neighbors
        cn4 = poit2[nind2].cnum;          //We get the cell to which it belongs

//We treat them as edges if the the neighbors are different
//And add the info on a site and cell level
//
	if ( (cn1 != cn2) && (cn3 != cn4) )
	{ //Hlay = Hlay + J12lay;
		getmatch = getmatch+1;	}


} // End the loop over neigbors



//	Hlay = Hlay + J1*nmatch;


	return getmatch;

}	//End the function


////////////////////////////////////////////////////////
//
////This will track the number of edges and introduce curvature

int bilayeredgecntK(Point poit[], Point poit2[], Cellprop CCC[], Cellprop PPP[],  int a, vector<int> n1, vector<int> n2)

{


int i, nind, nind2, getmatch = 0, cn1, cn2, cn3, cn4;

vector<int> layneigh;

//We only look at grid points who have neighbors in different cells
//

cn1 = poit[a].cnum;	//Cell number of first grid point

cn3 = poit2[a].cnum;             //Cell number of second grid point


//We loop over the hex neighbors
for (i=0; i<6; i++)

	{

	nind = n1[i];

	cn2 = poit[nind].cnum;

        nind2 = n2[i];        //We get the index of the neighbors
        cn4 = poit2[nind2].cnum;          //We get the cell to which it belongs

//We treat them as edges if the the neighbors are different
//And add the info on a site and cell level
//
	if ( (cn1 != cn2) && (cn3 != cn4) )
	{ getmatch = getmatch+1;	}


} // End the loop over neigbors

poit[a].nedge = getmatch;
poit[a].zhex = getmatch*CCC[cn1].kapp;

poit2[a].nedge = getmatch;
poit2[a].zhex = getmatch*PPP[cn2].kapp;



	return getmatch;

}	//End the function



///////////////////
//bookmark
//Here we keep track of the number of matching edges

double edgematchratavg(Cellprop CCC[], int a)
{

	int i =	1;

	int em = 0;

	for (i=1; i<a; i++)
	{em = em + CCC[i].edgemC;}

	if (a>0) { em = em/a;}


	return em;

}

////////////////////////////////////////////////////////
/////////////////////////////////////////////////
//Here we set up the inital neighbor list

void setncprops(Point poit[], Point poit2[],Cellprop CC[], Cellprop LL[], Grid& Ge, Grid& Le, Nlist& ne, Nlist& me, double L)
{

int tyl = Ge.Nnum;

int ss = Ge.NP;

int i = 0, ty1,ty2,e1,e2;

vector<int> c1, c2,b1,b2;

//Nlist li1, li2;

//li1.n = 0, li2.n = 0;

ne.List.clear();
me.List.clear();



int l1 = 0, l2 = 0;

int j = 0;

for (j=0; j<tyl; j++) {CC[j].CPclear(); LL[j].CPclear(); }

for (i=0; i<ss;i++)
{
	poit[i].Bi = 0;

	c1 = neighlistset(poit, Ge, CC,i);
	ty1 = poit[i].cnum; CC[ty1].Ai = CC[ty1].Ai + 1;
	CC[ty1].Pi = CC[ty1].Pi + c1[6];

	CC[ty1].Aihex = CC[ty1].Aihex + sqrt(3)*L*L/2;

	CC[ty1].Pihex = CC[ty1].Pihex + c1[6]*L/sqrt(3);

	//Here we check if it is a perimeter point and set up the neighbor list

	if (c1[6] > 0 ) {poit[i].Peri = 1; ne.List.push_back(i); l1 = l1+1;}
	else {poit[i].Peri = 0; }

	c2 = neighlistset(poit2, Le, LL,i);
        ty2 = poit2[i].cnum; LL[ty2].Ai = LL[ty2].Ai + 1;
        LL[ty2].Pi = LL[ty2].Pi + c2[6];

	LL[ty2].Aihex = LL[ty2].Aihex + sqrt(3)*L*L/2;

        LL[ty2].Pihex = LL[ty2].Pihex + c2[6]*L/sqrt(3);


        //Here we check if it is a perimeter point and set up the neighbor list

        if (c2[6] > 0 ) {poit2[i].Peri = 1; me.List.push_back(i); l2 = l2+1;}
	else {poit2[i].Peri = 0; }


	e1 =bilayeredgecntK(poit,poit2,CC,LL,i,c1,c2);
	e2 =bilayeredgecntK(poit2,poit,LL,CC,i,c2,c1);


	poit[i].nedge = e1;
	poit2[i].nedge = e2;

	CC[ty1].edgemC = CC[ty1].edgemC +e1;
	LL[ty2].edgemC = LL[ty2].edgemC +e2;

}

for (j=0; j<tyl; j++) 
{
	CC[j].Getmatchratio(); 
	LL[j].Getmatchratio(); 
}



ne.n = l1; me.n = l2;





}

///////////////////////////////////
///Here we calculate the geometric mean of activity of the cells for hex
//
//
double GMactivehex(Point poit[], int a)

{

	double GMa = poit[a].Activ;	//This will keep track of the GM
	int i = 0;
	int ind;
	int cel1, cel2;
	int acnt = 1;
	double aexp;
	double gma;
	cel1 = poit[a].cnum;

	for(i=0; i<6; i++)
	{ 
		ind = poit[a].neiind[i];
		cel2 = poit[ind].cnum;
	      	if (cel1 == cel2) {GMa = GMa*poit[ind].Activ; acnt=acnt+1; }
	}

	//For the geometric mean we take the nth root for n elements
	
	aexp = 1/acnt;

	gma = pow(GMa, aexp);

	return gma;
}





////////////////////////////////////////////////////////////
//bookmark
void getconnectivity( Point Po[], Grid G, Cellprop Cel[])
{


int nx = G.Lx;
int ny = G.Ly;
int nc = G.Nnum;
int La = G.Lat;

//The gap between the first and second cells is (nx+0.5)*Lat + 2*Lat


int i = 0, j = 0, k = 0;

int xcontin[nc][nx] = {0};
int ycontin[nc][ny] = {0};


//We'll initialize the discontinuity lines to a nonsense value like 

int x, y;
double Xhex, Yhex;

//First we do a left to right x direction sweep, going down columns (y direction ) then across

for (j=0; j<nx; j++)
{
        for (i=0; i<ny; i++)
        {
                for (k = 0; k< nc; k++ )
                {
                        if (Po[i*nx+j].cnum == k)
                                { xcontin[k][j] = 1; }  //Checking each column
                }       //Loop over cells

        } // Loop down the column

}       //Loop over the rows

//Then we do a y direction sweep, going down rows (x direction) then column


for (j=0; j<ny; j++)
{

        for (i=0; i<nx; i++)
        {
                for (k = 0; k< nc; k++ )
                {
                        if (Po[i+j*nx].cnum == k)
                                { ycontin[k][j] = 1; }  //Checking each row
                }       //Loop over cells

        } // Loop down the column

}       //Loop over the rows

//Now we have a record of whether each cell has a site in each column/row. We loop over the record and if there is a gap (ex. 1, 0, 1) we have a discontinuity


i = 0; j = 0;

for (j=0; j<nc; j++)
{
        Cel[j].connectx = 0; Cel[j].connecty = 0;
//First we check the boundaries
        if ((xcontin[j][0] == 0) &&(xcontin[j][1] ==1) &&(xcontin[j][nx-1] == 1) )
        {
                x = 1;
                Xhex = ((double)x+0.5)*La;
                //if (Laye == 2)
                //{ Xhex = ((double)x+0.5)*La +(nx+0.5)*La + 2*La; }
                Cel[j].connectx = 1; Cel[j].disconx = Xhex;
        }


if ((ycontin[j][0] == 0) &&(ycontin[j][1] ==1) &&(ycontin[j][ny-1] == 1) )
        {
                y = 1;
                Yhex = (double)y*La*sqrt(3)/2;
                Cel[j].connecty = 1; Cel[j].discony = Yhex;
        }

        if ((xcontin[j][nx-1] == 0) &&(xcontin[j][0] ==1) &&(xcontin[j][nx-2] == 1) )
        {
                x =(nx-1 % nx) + 1;
                Xhex = ((double)x+0.5)*La;
                //if (Laye == 2)
                //{Xhex = ((double)x+0.5)*La +(nx+0.5)*La + 2*La; }
                Cel[j].connectx = 1; Cel[j].disconx = Xhex;
        }
 if ((ycontin[j][ny-1] == 0) &&(ycontin[j][0] ==1) &&(ycontin[j][ny-2] == 1) )
        {
                y = floor((nx-1)/nx) + 1;
                Yhex = (double)y*La*sqrt(3)/2;
                Cel[j].connecty = 1; Cel[j].discony = Yhex;
        }
//Then the rest of the grid
        for(i=1; i<nx-1; i++)
        {
                if ((xcontin[j][i] == 0) &&(xcontin[j][i-1] ==1) &&(xcontin[j][i+1] == 1) )
                {
                x = (i % nx) + 1;
                Xhex = ((double)x+0.5)*La;
                //if (Laye == 2)
                //{ Xhex = ((double)x+0.5)*La +(nx+0.5)*La + 2*La; }
                Cel[j].connectx = 1; Cel[j].disconx = Xhex;

                }

	}

	for(i=1; i<ny-1; i++)
	{
       if ((ycontin[j][i] == 0) &&(ycontin[j][i-1] ==1) &&(ycontin[j][i+1] == 1) )
                {
                y = floor(i/nx) + 1;
                Yhex = (double)y*La*sqrt(3)/2;
                Cel[j].connecty = 1; Cel[j].discony = Yhex;
                }


        } //End the for loop over i



}       //End the for loop over j




}	//End the function
///////////////////////////////////////
//bookmark
//
//Here we get the x and y com of the cells
void getcom( Point Po[], Grid G, Cellprop Cel[], double Laye)
{

int nx = G.Lx;
int ny = G.Ly;
int le = G.NP;
int nc = G.Nnum;
int La = G.Lat;

double lenx = G.Lxhex;
double leny = G.Lyhex;

vector<double> xy = {0,0};

int i = 0;
int j = 0;
double xc[nc] = {0}, xcper[nc] = {0};
double yc[nc] = {0}, ycper[nc] = {0};
double xmin[nc] = {0}, xmax[nc] = {0};
double ymin[nc] = {0}, ymax[nc] = {0};
int npoint[nc] = {0};

double mdptx = G.COMx;
double mdpty = G.COMy;
int k = 0;

for (k = 0; k<nc; k++)
{ xmin[k] = 2*nx; xmax[k] = -2*nx; ymin[k] = 2*ny; ymax[k] = -2*ny; }

Point Go[le];

for (j=0; j<le; j++)
{
        Go[j] = Po[j];
        if (Laye == 2)
        {
                Go[j].xhex = Po[j].xhex -(nx+0.5)*La - 2*La;
        }

}


getconnectivity(Go,G,Cel);

//We need to account for the x values of the second layer

///////////////////////////////////////////////////////
//First we get the quadrant each com lies in

for (i=0; i<le; i++)
{
        for (j = 0; j< nc; j++ )
        {
                if (Go[i].cnum == j)
        { xc[j] = xc[j] + Go[i].xhex; yc[j] = yc[j] + Go[i].yhex;
                npoint[j] = npoint[j]+1;
//We also get the max and min of the x and y coordintates

                if (Go[i].xhex > xmax[j] ) { xmax[j] = Go[i].xhex; }

                if (Go[i].xhex < xmin[j] ) { xmin[j] = Go[i].xhex; }

                if (Go[i].yhex > ymax[j] ) { ymax[j] = Go[i].yhex; }

                if (Go[i].yhex < ymin[j] ) { ymin[j] = Go[i].yhex; }

        }       //End if for point

        }       //End loop over cells
}

///////////////////////////////////////////////////
//Now we work on the quadrant for each cell

i = 0; j = 0; k = 0;


for (k=0; k<nc; k++)
{
        if (npoint[k] > 0)
        {
                xc[k] =xc[k]/(double)npoint[k]; yc[k] = yc[k]/(double)npoint[k];
                if ( (xc[k] > mdptx ) && (yc[k] > mdpty ) ) { Cel[k].quad = 1;}

                if ( (xc[k] < mdptx ) && (yc[k] > mdpty ) ) { Cel[k].quad = 2;}

                if ( (xc[k] < mdptx ) && (yc[k] < mdpty ) ) { Cel[k].quad = 3;}

                if ( (xc[k] < mdptx ) && (yc[k] > mdpty ) ) { Cel[k].quad = 4;}


        } //Loop over midpoints
} //Loop over cells


////////////////////////////////////////////////////////////////
//Finally, we loop through the points one more time to get the adjusted coms

int lpnts = 0, lcels = 0;


for (lpnts=0; lpnts<le; lpnts++)
{
        for (lcels = 0; lcels< nc; lcels++ )
        {
                if (Go[lpnts].cnum == lcels)    {

        //First the simple case




                if (Cel[lcels].connectx == 0 ) {xcper[lcels] = xcper[lcels] + Go[lpnts].xhex;}
                if (Cel[lcels].connecty == 0 ) {ycper[lcels] = ycper[lcels] + Go[lpnts].yhex;}

        //Now the periodic cases        

                if (Cel[lcels].connectx == 1) {

//We adjust based on quadrants
                        if (Cel[lcels].quad == 1) {
                                if (Go[lpnts].xhex >= Cel[lcels].disconx ) {
                                        xcper[lcels] = xcper[lcels] + Go[lpnts].xhex; }
                                if (Go[lpnts].xhex < Cel[lcels].disconx) {
                                        xcper[lcels] = xcper[lcels] + Go[lpnts].xhex+lenx;}
                                        } //Quadrant 1
                        if (Cel[lcels].quad == 2) {
                                if (Go[lpnts].xhex > Cel[lcels].disconx ) {
                                        xcper[lcels] = xcper[lcels] +Go[lpnts].xhex-lenx; }
                                if (Go[lpnts].xhex <= Cel[lcels].disconx) {
                                        xcper[lcels] = xcper[lcels] + Go[lpnts].xhex;}
                                        } //Quadrant 2
                        if (Cel[lcels].quad == 3) {
                                if (Go[lpnts].xhex > Cel[lcels].disconx ) {
                                        xcper[lcels] = xcper[lcels] + Go[lpnts].xhex-lenx;}
                                if (Go[lpnts].xhex <= Cel[lcels].disconx) {
                                        xcper[lcels] = xcper[lcels] + Go[lpnts].xhex;}
                                        } //Quadrant 3
                        if (Cel[lcels].quad == 4) {
                                if (Go[lpnts].xhex >= Cel[lcels].disconx ) {
                                        xcper[lcels] = xcper[lcels] + Go[lpnts].xhex; }
                                if (Go[lpnts].xhex < Cel[lcels].disconx) {
                                        xcper[lcels] = xcper[lcels] + Go[lpnts].xhex+lenx;}
                                        } //Quadrant 4
                                } //X connectivity
                                ////////////////////////////////////////////////////////////////////    
                if (Cel[lcels].connecty == 1) {
//Now for the y 
                        if (Cel[lcels].quad == 1) {
                                if (Go[lpnts].yhex > Cel[lcels].discony ) {
                                        ycper[lcels] = ycper[lcels] + Go[lpnts].yhex-leny;}
                                if (Go[lpnts].yhex <= Cel[lcels].discony) {
                                        ycper[lcels] = ycper[lcels] + Go[lpnts].yhex;}
                                        } //Quadrant 1
                        if (Cel[lcels].quad == 2) {
                                if (Go[lpnts].yhex > Cel[lcels].discony ) {
                                        ycper[lcels] = ycper[lcels]+Go[lpnts].yhex-leny; }
                                if (Go[lpnts].yhex <= Cel[lcels].discony) {
                                        ycper[lcels] = ycper[lcels] + Go[lpnts].yhex;}
                                        } //Quadrant 2
                        if (Cel[lcels].quad == 3) {
                                if (Go[lpnts].yhex >= Cel[lcels].discony ) {
                                        ycper[lcels] = ycper[lcels] + Go[lpnts].yhex;}
                                if (Go[lpnts].yhex <= Cel[lcels].discony) {
                                        ycper[lcels] = ycper[lcels] + Go[lpnts].yhex+leny;}
                                        } //Quadrant 3
                        if (Cel[lcels].quad == 4) {
                                if (Go[lpnts].yhex >= Cel[lcels].discony ) {
                                        ycper[lcels] = ycper[lcels] + Go[lpnts].yhex; }
                                if (Go[lpnts].yhex < Cel[lcels].discony) {
                                        ycper[lcels] = ycper[lcels] + Go[lpnts].yhex+leny;}
                                        } //Quadrant 4

                                } //Y connectivity 

                }       //End if for point

        }       //End loop over cells
}       //End loop over points

//Then we divide by the number of elements
i = 0;
for (i=0; i<nc; i++)
{
        if (npoint[i] != 0 )
        {
                Cel[i].xcom = xcper[i]/(double)npoint[i];
                Cel[i].ycom = ycper[i]/(double)npoint[i];
		Cel[i].hexaorder = 0;
        }

        if (npoint[i] == 0 )
        {

                Cel[i].xcom = 0; Cel[i].ycom = 0;
		Cel[i].hexaorder = 0;
        }
}



}	//End the function

/////////////////////////////////////////////////////////
//This gets the COM at the very beginning of a simulation, assuming a confluent setup
//and none of the cells have crossed the periodic boundaries
/*
void getcominit( Point Po[], Grid G, Cellprop Cel[])
{

int nx = G.Lx;	int ny = G.Ly;	int le = G.NP;	int nc = G.Nnum; int La = G.Lat;

double lenx = G.Lxhex;
double leny = G.Lyhex;

int i = 0;
int j = 0;

int npoint[nc] = {0};

double mdptx = G.COMx,mdpty = G.COMy;

for (i=0; i<le; i++)
{



}

*/

////////////////////////////////////////////////

//Here we use the cc3d method to update the COM
//Xsite and Ysite will be the coordinates of the points, b will be 1 if adding and 0 if subtracting, and c will be the Cell area

void updatecomCC3D(Grid& G, Cellprop Cel[], double Xsite, double Ysite,int b,double c, int Cnum)
{

int sign;

double xglow = G.Xlo, xghi = G.Xhi, lenx = G.Lxhex, yglow = G.Ylo, yghi = G.Yhi, leny = G.Lyhex, gcomx = G.COMx, gcomy = G.COMy, xcomtemp = 0, ycomtemp = 0, La = G.Lat, rad3 = sqrt(3), Asite = La*La*rad3/2;

double xshift, xcellshift, pxshift, xnew, xchk, yshift, ycellshift,pyshift,ynew,ychk;

if (b == 1)
{
	sign = 1;
}

else
{
	sign = -1;
}



if (c != 0)
{

	xcomtemp = Cel[Cnum].xcom; ycomtemp = Cel[Cnum].ycom;
	xshift = xcomtemp - gcomx; yshift = ycomtemp - gcomy;
	xcellshift =xcomtemp -xshift; ycellshift =ycomtemp -yshift;
	pxshift = Xsite - xshift; pyshift = Ysite - yshift;
	xnew = (xcellshift*c+sign*pxshift)/(c+sign*Asite);
	ynew = (ycellshift*c+sign*pyshift)/(c+sign*Asite);
	xchk = xnew+xshift; ychk = ynew + yshift;
	xcomtemp = xchk;	ycomtemp = ychk;

	//Now we check if we are in the lattice

	if (xchk > xghi) 
	{

		xcomtemp = xchk - lenx;
	}


	if (xchk < xglow)
	{
		xcomtemp = xchk + lenx;
	}

	if (ychk > yghi)
	{
		ycomtemp = ychk - leny;
	}

	if (ychk < yglow)
	{
		ycomtemp = ychk + leny;
	}


}


else 
{
	if (b == 1)
	{
		xcomtemp = Xsite; ycomtemp = Ysite;
	}


	else
	{
		 Cel[Cnum].xcom = 0; Cel[Cnum].ycom = 0;

	}
}


	Cel[Cnum].xcom = xcomtemp; Cel[Cnum].ycom = ycomtemp;



}

/////////////////////////////////////////////////////////////////////
//This will get the com for all cells without changing the area

void getcomCC3D( Point Po[], Grid G, Cellprop Cel[])
{

int le = G.NP, nc = G.Nnum, i = 0, npoints[nc] = {0}, Cnum;

double xcomtemp = 0, ycomtemp = 0, La = G.Lat, rad3 = sqrt(3), Asite = La*La*rad3/2, Acel;


for (i=0; i<le; i++)
{
	Cnum = Po[i].cnum; xcomtemp = Po[i].xhex; ycomtemp = Po[i].yhex;
	Acel = npoints[Cnum]*Asite;
	updatecomCC3D(G,Cel,xcomtemp,ycomtemp,1,Acel, Cnum);
	npoints[Cnum]+=1;

}



}

///////////////////////////////////////////////////////////////////////
//Here we use a variation of the cc3d method to update the COM
//Xsite and Ysite will be the coordinates of the points, b will be 1 if adding and 0 if subtracting, and c will be the Cell area
//We look at the distance between the point and the cell COM, the point-Lx and the cell COM, and the point+Lx and the cell COM. We take the point with the minimum distance

void updatecommindist(Grid& G, Cellprop Cel[], double Xsite, double Ysite,int b,double c, int Cnum)
{

int sign;

double xglow = G.Xlo, xghi = G.Xhi, lenx = G.Lxhex, yglow = G.Ylo, yghi = G.Yhi, leny = G.Lyhex, gcomx = G.COMx, gcomy = G.COMy, La = G.Lat, rad3 = sqrt(3), Asite = La*La*rad3/2,xuse, yuse, xdist1, xdist2, xdist3, ydist1, ydist2, ydist3, xnew, ynew, cellcomx = Cel[Cnum].xcom, cellcomy = Cel[Cnum].ycom;


if (b == 1)
{
	sign = 1;
}

else
{
	sign = -1;
}

//First the X coordinate

xdist1 = (cellcomx-Xsite)*(cellcomx-Xsite);
xuse = Xsite;

xdist2 = (cellcomx-(Xsite-lenx))*(cellcomx-(Xsite-lenx));
if (xdist2 < xdist1) 
{
	xuse = Xsite-lenx; 
}

xdist3 = (cellcomx-(Xsite+lenx))*(cellcomx-(Xsite+lenx));

if ( (xdist3 < xdist2) && (xdist3 < xdist1) )
{
	xuse = Xsite + lenx;
}


ydist1 = (cellcomy-Ysite)*(cellcomy-Ysite);
yuse = Ysite;

ydist2 = (cellcomy-(Ysite-leny))*(cellcomy-(Ysite-leny));
if (ydist2 < ydist1)
{
        yuse = Ysite-leny; 
}

ydist3 = (cellcomy-(Ysite+leny))*(cellcomy-(Ysite+leny));

if ( (ydist3 < ydist2) && (ydist3 < ydist1) )
{
        yuse = Ysite + leny;
}


if (c != 0)
{

	xnew = (cellcomx*c+sign*xuse*Asite)/(c+sign*Asite);
	ynew = (cellcomy*c+sign*yuse*Asite)/(c+sign*Asite);

	//Now we check if we are in the lattice


}


else 
{
	if (b == 1)
	{
		xnew = Xsite; ynew = Ysite;
	}


	else
	{
		 Cel[Cnum].xcom = 0; Cel[Cnum].ycom = 0;

	}
}


	Cel[Cnum].xcom = xnew; Cel[Cnum].ycom = ynew;



}

/////////////////////////////////////////////////////////////////////
//This will get the com for all cells without changing the area

void getcommindist( Point Po[], Grid G, Cellprop Cel[])
{

int le = G.NP, nc = G.Nnum, i = 0, npoints[nc] = {0}, Cnum;

double xcomtemp = 0, ycomtemp = 0, La = G.Lat, rad3 = sqrt(3), Asite = La*La*rad3/2, Acel;


for (i=0; i<nc; i++)
{
	Cel[i].xcom = 0; Cel[i].ycom = 0;
}


for (i=0; i<le; i++)
{
	Cnum = Po[i].cnum; xcomtemp = Po[i].xhex; ycomtemp = Po[i].yhex;
	Acel = npoints[Cnum]*Asite;
	updatecommindist(G,Cel,xcomtemp,ycomtemp,1,Acel, Cnum);
	npoints[Cnum]+=1;

}



}





////////////////////////////////////////////////////////////////
//Now we start looking at hexatic order
//
//Here we get the x and y com of the cells
void hexorder( Point Po[], Cellprop CC[], Grid G )
{

int le = G.NP;
int ncel = G.Nnum;


int i = 0,k = 0;
double rx = 0,ry = 0;
double theta = 0;
int nc[ncel] = {0};
std::complex<double> phase;
std::complex<double> psi[ncel];

std::vector<Cellprop> LL;
LL.resize(ncel);

for (i=0;i<ncel; i++)
{	psi[i] = std::complex<double>(0,0);  //LL[i] = CC[i];
//	if (laye == 2) 
//		{LL[i].xcom = CC[i].xcom + (nx+0.5)*La + 2*La; } 
}

for (i=0; i<le; i++)
{
	for (k=0; k<ncel; k++)	
		{
			if (Po[i].cnum == k) 
				{ rx = Po[i].xhex-CC[k].xcom; ry = Po[i].yhex-CC[k].ycom; nc[k] = nc[k]+1; theta = atan2(ry,rx); phase = std::complex<double>(0.0,6*theta); psi[k] += std::exp(phase);   }	
		}

}

//Then we average

for (k=0; k<ncel; k++) 
{
	if (nc[k] > 0)  {CC[k].hexaorder = psi[k].real() / nc[k];} 
	if (nc[k] == 0) { CC[k].hexaorder = 0; }

} //End loop over cells


}

//////////////////////////////////////////////////
/*

void hexorder2( Cellprop CC[], Grid G, int a )
{

int le = G.NP;
int ncel = G.Nnum;


int i = 0,k = 0;
double rx = 0,ry = 0;
double theta = 0;
int nc[ncel] = {0};
std::complex<double> phase;
std::complex<double> psi[ncel];

std::vector<Cellprop> LL;
LL.resize(ncel);

for (i=0;i<ncel; i++)
{	psi[i] = std::complex<double>(0,0); }

for (i=0; i<le; i++)
{
	for (k=0; k<ncel; k++)	
		{
			if (Po[i].cnum == k) 
				{ rx = Po[i].xhex-CC[k].xcom; ry = Po[i].yhex-CC[k].ycom; nc[k] = nc[k]+1; theta = atan2(ry,rx); phase = std::complex<double>(0.0,6*theta); psi[k] += std::exp(phase);   }	
		}

}

//Then we average

for (k=0; k<ncel; k++) 
{
	if (nc[k] > 0)  {CC[k].hexaorder = psi[k].real() / nc[k];} 
	if (nc[k] == 0) { CC[k].hexaorder = 0; }

} //End loop over cells


}


*/

////////////////////////////////////
//
double meanhexorder( Cellprop CC[], Grid G)
{


int ncel = G.Nnum;

int i = 0;

double meanpsi = 0; 

for (i=0; i<ncel; i++)
{

	meanpsi += CC[i].hexaorder;

}


if (ncel > 0) { meanpsi /= ncel; }

return meanpsi;

}



///////////////////////////////////
//Here we get the com in hex coordinates
//

double getadhinlay(Point poit[], Cellprop CC[], int a)
{


	double Had = 0;
	int cenum = poit[a].cnum;
	double J1 = CC[cenum].Ji;
	double J12;
	int i;
	int idn,cn;

	for (i=0; i<6; i++)

	{
		idn = poit[a].neiind[i];
		cn = poit[idn].cnum;

		if (cn != cenum) { J12 = (J1 + CC[cn].Ji)/2; Had+=J12; }

	}
		

	return Had;
}

///////////////////////////////////////////////////
//This will setup the inital grid


void setgrid(Point poit[], Grid gee, Cellprop CC[], int lay )
{

int nps = gee.NP;
int nall = gee.Nnum;
std::uniform_int_distribution<uint32_t> uint_distcell(0,nall);


int rnum;
int i = 0;
int j = 0;
int ty = 0;

for (j=0; j<nall; j++ ) {CC[j].CPclear();}

for (i = 0; i < nps; i++ )

{



        poit[i].initpts(lay);

        poit[i].cnum = 0;

	poit[j].Peri = 0;

        rnum = uint_distcell(rng);

        ty = 1;

        for (ty = 1; ty < nall; ty++ )

        {
                if ( (CC[ty].Ai < CC[ty].Ainit ) && (rnum == ty)  )

                { poit[i].cnum = ty;  CC[ty].Ai = CC[ty].Ai + 1;}

        }



}       //ends the loop over the grid

//Now we need to make sure we didn't underassign


int tchk = 0;
int ty2 = 0;
int vc;

for (j = 0; j < nps; j++ )

{

        ty2 = 1;

        vc = poit[j].cnum;

        tchk = 0;

        for (ty2 = 1; ty2 < nall; ty2++ )

        {

               if ( (CC[ty2].Ai < CC[ty2].Ainit ) && (vc == 0) && (tchk == 0 ))

                { poit[j].cnum = ty2; CC[ty2].Ai = CC[ty2].Ai+1; tchk = tchk+1;}

        }


} //End the second loop


}


/////////////////////////////////////////
//
//This will setup the inital grid with a less random setup


void setgrid2(Point poit[], Grid gee, Cellprop CC[], int lay )
{

int nps = gee.NP;
int nall = gee.Nnum;
std::uniform_int_distribution<uint32_t> uint_distcell(0,nall);


int i = 0;
int j = 0;
int num = 0;


for (j=0; j<nall; j++ ) {CC[j].CPclear();}

for (i = 0; i < nps; i++ )

{

        poit[i].initpts(lay);
 
	poit[i].Peri = 0;

	if (CC[num].Ai<= CC[num].Ainit )
	{ poit[i].cnum = num; CC[num].Ai = CC[num].Ai + 1; }

	if (CC[num].Ai == CC[num].Ainit )
	{ num = num+1; }

	if ( (CC[nall-1].Ai == CC[nall-1].Ainit ) && ( i < nps ) )

	{ poit[i].cnum = 0; CC[0].Ai = CC[0].Ai + 1; }
	
}       //ends the loop over the grid



}	//Ends the function

/////////////////////////////////////
//We do a similar non random setup for a confluent system
//

void setgrid3(Point poit[], Grid gee, vector<int> nc, int lay )
{

int nps = gee.NP;
int ncells = gee.Nnum-1;


int i = 0,num = 0, cel = 0;


for (i = 0; i < nps; i++ )

{

        poit[i].layer = lay;

        poit[i].Peri =  0;

        poit[i].Activ = 0;
 

	if (cel == ncells)
	{
		poit[i].cnum = ncells;
	}


	if (cel < ncells )
	{
		
		if (num == nc[cel])
		{ cel = cel +1; num = 0; }

		if (num < nc[cel] )
		{ num = num+1; poit[i].cnum = cel; }


	}


}       //ends the loop over the grid



}	//Ends the function

///////////////////////////////////////////////////////
//This will set up the grid in the x direction for one layer and the y direction for the other

void setgrid4(Point poit[], Grid gee, vector<int> nc, int lay )
{

int nps = gee.NP;
int ncells = gee.Nnum;
int nx = gee.Lx;
int ny = gee.Ly;


int i = 0,j = 0, index = 0,num = 0, cel = 0, extra = nps;

std::uniform_int_distribution<uint32_t> cellextra(2,ncells-2);

int rngcell = cellextra(rng);

int nc2[ncells+1];

for (i=0; i<ncells; i++)
{nc2[i] = nc[i]; extra -= nc[i]; }

nc2[rngcell] = nc2[rngcell]+ extra;


for (i = 0; i < nx; i++ )	//Loop down rows
{

	for (j=0; j<ny; j++)	//Loop down columns 
	{ 

		//if (lay == 2)
		//	{index = i + j*nx; }
		//if (lay != 2)
		//	{index = i*nx+j; }

		index = i + j*nx;
        	poit[index].layer = lay;
        	poit[index].Peri =  0;
        	poit[index].Activ = 0;

		if (cel < ncells )
		{
		
			if (num == nc2[cel])
			{ cel = cel +1; num = 0; }

			if ( (num < nc2[cel] ) && ( cel < ncells ) )
			{ num = num+1; poit[index].cnum = cel; }


		}


	}       //ends the loop down columns



}	//Ends the loop down rows



}	//Ends the function




///////////////////////////////////////////////////////////////
//Here we set up the list of moves that could change cell

void getnlist(Point poit[], Grid gee, Nlist nei)
{
int nps = gee.NP;

int i = 0;
int j = 0;


nei.List.clear();

for (i = 0; i<nps; i++)
{ if (poit[i].Peri == 1) { j = j+1; nei.List.push_back(i); }	}

nei.n = j;

}

//////////////////////////
//Here we setup the grid using a discrete Voronoi algorithm

void setgridVoro(Point pt[], Grid Ge, int empty, int lay)

{
int nps = Ge.NP;
int nall = Ge.Nnum;

int Lenx = Ge.Lx;
int Leny = Ge.Ly;

//std::uniform_int_distribution<uint32_t> uint_initx(0,100*Lenx);

//std::uniform_int_distribution<uint32_t> uint_inity(0,100*Leny);

std::uniform_int_distribution<uint32_t> uint_shift(0,360);

std::uniform_int_distribution<uint32_t> uint_scale(1,5);

std::uniform_int_distribution<uint32_t> uint_distcell(1,nall-1);



int rngx, rngy, rngang, rnscl;

//We get the com of the grid
//
double xgcom = Lenx/2, ygcom = Leny/2;


//double crad = Lenx/5;		//Radius of the unit sphere

double crad = (Lenx+Leny)/6;

//If we want N cells we pick N random points

int i = 0;

int j = 0;


double dchkmin;

double dist; 

double xcel[nall], ycel[nall];

double rshift;

rngang = uint_shift(rng);

rshift = (double)rngang*M_PI/180;


for (i = 0; i < nall; i++)
{

//        rngx = uint_initx(rng);

//        rngy = uint_inity(rng);

	rnscl = uint_scale(rng);


//We use the points to parametrize a sphere near the center
	

	xcel[i] = xgcom + crad*cos(2*M_PI*i/nall + rshift)/((double)rnscl);

	ycel[i] = ygcom + crad*sin(2*M_PI*i/nall + rshift)/((double)rnscl);	
	
	

}



//Now we loop through the sites and organize them based on distance

int k = 0;

int site = 0;

for (j=0; j<nps; j++)
{
	dchkmin = (Lenx+1)*(Leny+1);
	for (k=0; k<nall; k++)
	{
		dist = (pt[j].xhex-xcel[k])*(pt[j].xhex-xcel[k])+
		       (pt[j].yhex-ycel[k])*(pt[j].yhex-ycel[k]);

		dist = sqrt(dist);

		if (dist < dchkmin) {dchkmin = dist; site = k; }

	}

	//The site with the minimum distance is the cell
	pt[j].cnum = site;	

	if (empty == 1 )

	{ 
		if ( site==0 ) {pt[j].cnum =1;}

//		if ( site==extra-1 ) {pt[j].cnum = extra; }

	}


	pt[j].layer = lay;

	pt[j].Peri =  0;

	pt[j].Activ = 0;
}



//end the function
//
}

///////////////////////////////////////////////////////
//Here we set up the grid by using random points on the lattice as Voronoi seeds
void setgridVoro2(Point pt[], Grid Ge, int empty, int lay)
{

int nps = Ge.NP;
int nall = Ge.Nnum;

int Lenx = Ge.Lx;
int Leny = Ge.Ly;

//std::uniform_int_distribution<uint32_t> uint_initx(0,100*Lenx);

//std::uniform_int_distribution<uint32_t> uint_inity(0,100*Leny);

//std::uniform_int_distribution<uint32_t> uint_shift(0,360);

//std::uniform_int_distribution<uint32_t> uint_scale(1,5);

//std::uniform_int_distribution<uint32_t> uint_distcell(1,nall-1);

std::uniform_int_distribution<uint32_t> uint_distpoints(0,nps);

vector<int> rsites;

int i, j, rindx, nsites = 0, nck;

double dchkmin;

double dist; 

double xcel[nall], ycel[nall];


for (i = 0; i < nall; i++)
{

	rindx = uint_distpoints(rng);

	if (nsites == 0) 
		{ 
			rsites.push_back(rindx); 
	
			xcel[i] = pt[rindx].xhex;

        		ycel[i] = pt[rindx].yhex;

			nsites = nsites+1;
	
		}


	if (nsites >0 )
	{

		nck = 0;
		for (j=0; j<nsites; j++)
		{
			if (rsites[j] == rindx) 
			{
				nck = 1; i-=1;
			}

		}

		if (nck == 0)
		{

			rsites.push_back(rindx);

                        xcel[i] = pt[rindx].xhex;

                        ycel[i] = pt[rindx].yhex;

                        nsites = nsites+1;

		}


	}	

}



//Now we loop through the sites and organize them based on distance

int k = 0;

int site = 0;

for (j=0; j<nps; j++)
{
	dchkmin = (Lenx+1)*(Leny+1);
	for (k=0; k<nall; k++)
	{
		dist = (pt[j].xhex-xcel[k])*(pt[j].xhex-xcel[k])+
		       (pt[j].yhex-ycel[k])*(pt[j].yhex-ycel[k]);

		dist = sqrt(dist);

		if (dist < dchkmin) {dchkmin = dist; site = k; }

	}

	//The site with the minimum distance is the cell
	pt[j].cnum = site;	

	if (empty == 1 )

	{ 
		if ( site==0 ) {pt[j].cnum =1;}

//		if ( site==extra-1 ) {pt[j].cnum = extra; }

	}


	pt[j].layer = lay;

	pt[j].Peri =  0;

	pt[j].Activ = 0;
}





}

//////////////////////////////////////////////
//bookmark
//Here we run the sorting algorithm using a hex lattice

void CPsorthex(Point poit[], Grid gee, Cellprop CC[], Hamil engy, double Temp)  

{								

//We get the number of points in the system

int nps = gee.NP;


std::uniform_int_distribution<uint32_t> uint_distAll(0,nps-1);

std::uniform_int_distribution<uint32_t> uint_prob(1,100);



engy.Hclear();

//We randomly pick a point on the grid

int ranum;

ranum = uint_distAll(rng);

//We get the cell the point belongs to

int idc1;

idc1 = poit[ranum].cnum;

//We pick a different point


int ranum2;

ranum2 = uint_distAll(rng);

int idc2;

//We get the cell this belongs to

idc2 = poit[ranum2].cnum; 

//We create a point array for the switch

Point ptmv[nps];

int i = 0;

for (i=0; i<nps; i++) {ptmv[i] = poit[i]; }

//Then we switch the two chosen points in the new structure

ptmv[ranum2] = poit[ranum];

ptmv[ranum] = poit[ranum2];


int hrow1 = poit[ranum].hexrw;
int hrow2 = poit[ranum2].hexrw;




//If the points belong to different cells we switch		
		
if  ( idc1 != idc2  )
{

//We make a vector to get the neighbors
//
vector<int> neigh;

//neigh = neighptslistMoor(poit, gee, CC,ranum);
neigh = neighlistset(poit,gee,CC, ranum);


//We calculate the energy due to adhesion
//
double HH;


HH = getadhinlay(poit, CC, ranum);

engy.Hadhold = engy.Hadhold + HH;


////////////////////

vector<int> neigh2;

//neigh2 = neighptslistMoor(poit, gee,CC,ranum2);
neigh2 = neighlistset(poit,gee,CC, ranum2);

HH = getadhinlay(poit, CC, ranum2);

engy.Hadhold = engy.Hadhold + HH;


//Now we move the grid points
//////////////////////
//bookmark switch



neigh.clear();
//neigh = neighptslistMoor(ptmv, gee, CC,ranum);
neigh = neighlistset(ptmv,gee,CC, ranum);


HH = getadhinlay(ptmv, CC, ranum);

engy.Hadhnew = engy.Hadhnew + HH;

////////////////////

neigh2.clear();
//neigh2 = neighptslistMoor(ptmv, gee, CC,ranum2);
neigh2 = neighlistset(ptmv,gee,CC, ranum2);

HH = getadhinlay(ptmv, CC, ranum2);

engy.Hadhnew = engy.Hadhnew + HH;

////////////////

//Now we do the energy check
//


double Hnew = engy.Hadhnew;
double Hold = engy.Hadhold;
double Bprob, msprob;

Bprob = exp( -(Hnew-Hold)/Temp );



// And now we check
//
		
if (Hnew < Hold )
{poit[ranum] = ptmv[ranum]; poit[ranum2] = ptmv[ranum2];
poit[ranum].hexrw = hrow2; poit[ranum2].hexrw = hrow1; }	

	
if (Hnew >= Hold )

{

	msprob = (double)uint_prob(rng)/100;

	if ( Bprob >= msprob)

		{poit[ranum] = ptmv[ranum]; poit[ranum2] = ptmv[ranum2];
		poit[ranum].hexrw = hrow2; poit[ranum2].hexrw = hrow1;}


}

		

//ptmv[ranum] = poit[ranum];
//ptmv[ranum2] = poit[ranum2];

} // This ends the case of two different cell types



}

//////////////////////////////////////////////////////////

//bookmark
//Here we run an optimized sorting algorithm

void CPsorthex2(Point poit[], Grid gee, Cellprop CC[], Hamil engy, double Temp, Nlist& nl)  
{								

//We get the number of points in the system

int nps = gee.NP;


int Nperi = nl.n;

std::uniform_int_distribution<uint32_t> uint_distAll(0,nps-1);

std::uniform_int_distribution<uint32_t> uint_Nperi(0,Nperi-1);

std::uniform_int_distribution<uint32_t> uint_prob(1,100);



engy.Hclear();

//We randomly pick a point on the grid that can change

int ranum0;

int ranum;

ranum0 = uint_Nperi(rng);

//ranum = uint_distAll(rng);

ranum = nl.List[ranum0];

//We get the cell the point belongs to

int idc1;

idc1 = poit[ranum].cnum;

//We pick a different point


int ranum1;

ranum1 = uint_Nperi(rng);


int ranum2;

//ranum2 = uint_distAll(rng);

ranum2 = nl.List[ranum1];

int idc2;

//We get the cell this belongs to

idc2 = poit[ranum2].cnum; 

//We create a point array for the switch

Point ptmv[nps];

int i = 0;

for (i=0; i<nps; i++) {ptmv[i] = poit[i]; }




//Then we switch the two chosen points in the new structure

ptmv[ranum2].cnum = poit[ranum].cnum;

ptmv[ranum].cnum = poit[ranum2].cnum;

//We redo if the cells are the same



//If the points belong to different cells we switch		
		
if  ( idc1 != idc2  )
{

//We make a vector to get the neighbors
//
vector<int> neigh;
vector<int> neigh1neighs;

int ic1;

//neigh = neighptslistMoor(poit, gee, CC,ranum);
neigh = neighlistset(poit,gee,CC, ranum);

if (neigh[6] > 0 ) {poit[ranum].Peri = 1;}
if (neigh[6] == 0 ) {poit[ranum].Peri = 0;}

i = 0; for (i = 0; i<6; i++) {ic1 = neigh[i]; neigh1neighs.push_back(ic1);}



//We calculate the energy due to adhesion
//
double HH;


HH = getadhinlay(poit, CC, ranum);

engy.Hadhold = engy.Hadhold + HH;


////////////////////

vector<int> neigh2;
vector<int> neigh2neighs;
int ic2;

//neigh2 = neighptslistMoor(poit, gee,CC,ranum2);
neigh2 = neighlistset(poit,gee,CC, ranum2);

if (neigh2[6] > 0 ) {poit[ranum2].Peri = 1;}
if (neigh2[6] == 0 ) {poit[ranum2].Peri = 0;}


i = 0; for (i = 0; i<6; i++) {ic2 = neigh2[i]; neigh2neighs.push_back(ic2);}



HH = getadhinlay(poit, CC, ranum2);

engy.Hadhold = engy.Hadhold + HH;


//Now we move the grid points
//////////////////////
//bookmark switch



neigh.clear();
//neigh = neighptslistMoor(ptmv, gee, CC,ranum);
neigh = neighlistset(ptmv,gee,CC, ranum);

if (neigh[6] > 0 ) {ptmv[ranum].Peri = 1;}
if (neigh[6] == 0 ) {ptmv[ranum].Peri = 0;}


HH = getadhinlay(ptmv, CC, ranum);

engy.Hadhnew = engy.Hadhnew + HH;

////////////////////

neigh2.clear();
//neigh2 = neighptslistMoor(ptmv, gee, CC,ranum2);
neigh2 = neighlistset(ptmv,gee,CC, ranum2);

if (neigh2[6] > 0 ) {ptmv[ranum2].Peri = 1;}
if (neigh2[6] == 0 ) {ptmv[ranum2].Peri = 0;}

HH = getadhinlay(ptmv, CC, ranum2);

engy.Hadhnew = engy.Hadhnew + HH;

////////////////

//Now we do the energy check
//


double Hnew = engy.Hadhnew;
double Hold = engy.Hadhold;
double Bprob, msprob;

Bprob = exp( -(Hnew-Hold)/Temp );

vector<int>nneighs1, nneighs2;

// And now we check
//
		
if (Hnew < Hold )
{poit[ranum].cnum = ptmv[ranum].cnum; poit[ranum2].cnum = ptmv[ranum2].cnum;
poit[ranum].Peri = ptmv[ranum].Peri; poit[ranum2].Peri = ptmv[ranum2].Peri; 

//auto k1;

//auto k2;

i = 0; for (i = 0; i<6; i++) {
	
ic1 = neigh1neighs[i];	ic2 = neigh2neighs[i];

nneighs1 = neighlistset(ptmv,gee,CC, ic1);

if (nneighs1[6] > 0) { poit[ic1].Peri = 1; 

	auto k1 = find(nl.List.begin(), nl.List.end(),ic1);

        if ( k1!= nl.List.end()){nl.List.erase(nl.List.begin(),k1); nl.n = nl.n-1;}

         nl.List.push_back(ic1); nl.n = nl.n+1;

	}

nneighs2 = neighlistset(ptmv,gee,CC, ic2);

if (nneighs2[6] > 0 )
	{
	poit[ic2].Peri = 1;

        auto k2 = find(nl.List.begin(), nl.List.end(),ic2);

        if (k2!=nl.List.end()){nl.List.erase(nl.List.begin(),k2); nl.n = nl.n-1;}

         nl.List.push_back(ic2); nl.n = nl.n+1;


	}



	}	//ends the neighbor loop


}	//ends the case of lower energy

//poit[ranum].hexrw = hrow2; poit[ranum2].hexrw = hrow1; }	

	
/////////////////////

if (Hnew >= Hold )

{

	msprob = (double)uint_prob(rng)/100;

	if ( Bprob >= msprob)

{poit[ranum].cnum = ptmv[ranum].cnum; poit[ranum2].cnum = ptmv[ranum2].cnum; 
poit[ranum].Peri = ptmv[ranum].Peri; poit[ranum2].Peri = ptmv[ranum2].Peri;

//auto k1;

//auto k2;

i = 0; for (i = 0; i<6; i++) {
	
ic1 = neigh1neighs[i];	ic2 = neigh2neighs[i];

nneighs1 = neighlistset(ptmv,gee,CC, ic1);

if (nneighs1[6] > 0) { poit[ic1].Peri = 1; 

	auto k1 = find(nl.List.begin(), nl.List.end(),ic1);

        if (k1<nl.List.end()){nl.List.erase(nl.List.begin(),k1); nl.n = nl.n-1;}

         nl.List.push_back(ic1); nl.n = nl.n+1;

	}

nneighs2 = neighlistset(ptmv,gee,CC, ic2);

if (nneighs2[6] > 0 )
	{
	poit[ic2].Peri = 1;

        auto k2 = find(nl.List.begin(), nl.List.end(),ic2);

        if (k2<nl.List.end()){nl.List.erase(nl.List.begin(),k2); nl.n = nl.n-1;}

         nl.List.push_back(ic2); nl.n = nl.n+1;

	}




	} //ends the neighbor loop
//		poit[ranum].hexrw = hrow2; poit[ranum2].hexrw = hrow1;}


}	//ends the successful prob
	//
	//

} 	//ends the prob check

	


//ptmv[ranum] = poit[ranum];
//ptmv[ranum2] = poit[ranum2];

} // This ends the case of two different cell types


}



//////////////////////////////////////


void setbendingy(Cellprop Ce[], Bendlist &blist, int ncell, double reg1, double reg2, double Bip, double Kengy)
{


	int celloop = 1;

	blist.nreg1 = 0; blist.nreg2 = 0;

	for (celloop = 1; celloop < ncell; celloop++)
	{
		if (Ce[celloop].ycom <= reg1 )
		{
			Ce[celloop].Jlay = Bip; 
			blist.reg1list.push_back(celloop);
			blist.nreg1+=1;
		}

		else if ((Ce[celloop].ycom > reg1) && (Ce[celloop].ycom <= reg2))	
		{
			Ce[celloop].Jlay = Bip - Kengy;
			blist.reg2list.push_back(celloop);
			blist.nreg2+=1;
		}



	}



}

////////////////////////////////////////////////////////////////////
void setbendingx(Cellprop Ce[], Bendlist &blist, int ncell, double reg1, double reg2, double Bip)
{


	int celloop = 1;

	blist.nreg1 = 0; blist.nreg2 = 0;

	for (celloop = 1; celloop < ncell; celloop++)
	{
		if (Ce[celloop].xcom <= reg1 )
		{
			Ce[celloop].Jlay = Bip; 
			blist.reg1list.push_back(celloop);
			blist.nreg1+=1;
		}

		else if ((Ce[celloop].xcom > reg1) && (Ce[celloop].xcom <= reg2))	
		{
			Ce[celloop].Jlay = Bip/2;
			blist.reg2list.push_back(celloop);
			blist.nreg2+=1;
		}



	}



}

/////////////////////////////////////////////////////
//This will return the number of cells with high matching ratio in the region
double regionbimatchy(Cellprop Ce[],int ncell, double reg1, double reg2)
{

	int cel=1, nmatch = 0;
	double coord, matchavg = 0;

	for (cel=1; cel<ncell; cel++)
	{
		coord = Ce[cel].ycom;

		if ( ( coord < reg2) && ( coord >= reg1 ) ) 
			{	
				nmatch +=1;
				matchavg += Ce[cel].matchratio;
			}

	}

	if (nmatch != 0)
		{matchavg = matchavg/(double)nmatch; }
	else
		{std::cout << "There were " << nmatch << " cells in the region so we move," << std::endl;
		matchavg = 1;}

	return matchavg;


}


////////////////////////////////////////////////////////////
double regionbimatchx(Cellprop Ce[],int nal, double reg1)
{

	int cel=1, nmatch = 0;
	double coord, matchavg;

	for (cel=1; cel<nal; cel++)
	{
		coord = Ce[cel].xcom;

		if ( coord <= reg1)
			{	
				nmatch +=1;
			 	matchavg += Ce[cel].matchratio;
			}

	}


	if (nmatch != 0)
		{matchavg = matchavg/(double)nmatch; }
	else
		{matchavg = 1; }

	return matchavg;


}


//////////////////////////////////////////////////
//bookmark

void setbendingpointy(Point Po[], int npts, double y1, double y2, double y3, double y4, double Bip, double Kengy)
{


        int poitloop=0,newpointslow = 0, newpointsmid = 0, newpointshi = 0;
        double coord;
	int c1 = 0, c2 = 0;

	for (poitloop = 0; poitloop < npts; poitloop++)
	{
		coord = Po[poitloop].yhex;


		if ( (coord >= y1 ) && (coord < y2 ) )
                {

                        Po[poitloop].Bi = Bip - Kengy; newpointslow+=1;
                }


		else if ( (coord >= y2 ) && (coord < y3) )
		{
			
			Po[poitloop].Bi = Bip; newpointsmid+=1;
			
		}

		else if ((coord >= y3 ) && (coord < y4))	
		{
			Po[poitloop].Bi = Bip - Kengy; newpointshi+=1;
		}

	}


	std::cout << "There are " << newpointslow << " points in the lower region, " << newpointsmid << " points in the middle region, and  " << newpointshi << " points in the top region. " << std::endl;


}


////////////////////////////////////////////////////////////
//

void setbendingpointysharp(Point Po[], int npts, double y1, double y2, double Bip)
{


        int poitloop=0,newpointsmid = 0;
        double coord;

	for (poitloop = 0; poitloop < npts; poitloop++)
	{
		coord = Po[poitloop].yhex;


		if ( (coord >= y1 ) && (coord < y2) )
		{
			
			Po[poitloop].Bi = Bip; newpointsmid+=1;
			
		}


	}




}




////////////////////////////////////////////////////////////////////
void setbendingpointx(Point Po[], int npts, double x1, double x2, double x3, double x4, double Bip, double Kengy)
{

   	int poitloop=0,newpointslow = 0, newpointsmid = 0, newpointshi = 0;
        double coord;
	int c1 = 0, c2 = 0;

	for (poitloop = 0; poitloop < npts; poitloop++)
	{
		coord = Po[poitloop].xhex;


		if ( (coord >= x1 ) && (coord < x2 ) )
                {

                        Po[poitloop].Bi = Bip - Kengy; newpointslow+=1;
                }


		else if ( (coord >= x2 ) && (coord < x3) )
		{
			
			Po[poitloop].Bi = Bip; newpointsmid+=1;
			
		}

		else if ((coord >= x3 ) && (coord < x4))	
		{
			Po[poitloop].Bi = Bip - Kengy; newpointshi+=1;
		}

	}


	std::cout << "There are " << newpointslow << " points in the lower region, " << newpointsmid << " points in the middle region, and  " << newpointshi << " points in the top region. " << std::endl;




}


//////////////////////////////////////////////////
//This will introduce bending in a more straightfoward way
void siteneighheight(Point Po[], int a )
{

	int neiind, neighi, indloop = 0;
	double avgheight = 0;

	for (indloop=0; indloop<6; indloop++)
	{
		neiind = Po[a].neiind[indloop];
		avgheight = Po[neiind].nedge;
	}



}



/////////////////////////////////////////////////////
//This will return the number of cells with high matching ratio in the region
double regionbimatchpointy(Point Po[], Nlist& nei,double reg1, double reg2)
{

	int poit=0, nmatch = 0, npoints = nei.n, index;
	double coord, matchavg = 0;

	for (poit=0; poit<npoints; poit++)
	{
		index = nei.List[poit];
		coord = Po[index].yhex;

		if ( ( coord < reg2) && ( coord >= reg1 ) ) 
			{	
				nmatch +=1;
				matchavg += (double)Po[index].nedge/(double)Po[index].nperi;
			}

	}

	if (nmatch != 0)
		{matchavg = matchavg/((double)nmatch); }
	else
		{std::cout << "There were " << nmatch << " points in the region so we move," << std::endl;
		matchavg = 1;}

	return matchavg;


}

/////////////////////////////////////////////////
//This will include both endpoints
double regionbimatchpointy2(Point Po[], Nlist& nei,double reg1, double reg2)
{

	int poit=0, nmatch = 0, npoints = nei.n, index;
	double coord, matchavg = 0;

	for (poit=0; poit<npoints; poit++)
	{
		index = nei.List[poit];
		coord = Po[index].yhex;

		if ( ( coord <= reg2) && ( coord >= reg1 ) ) 
			{	
				nmatch +=1;
				matchavg += (double)Po[index].nedge/(double)Po[index].nperi;
			}

	}

	if (nmatch != 0)
		{matchavg = matchavg/((double)nmatch); }
	else
		{std::cout << "There were " << nmatch << " points in the region so we move," << std::endl;
		matchavg = 1;}

	return matchavg;


}




////////////////////////////////////////////////////////////////
double newregionbimatchpointy(Point Po[], Nlist& nei,double reg1, double reg2, double reg3, double reg4)
{

	int poit=0, nmatch = 0, npoints = nei.n, index;
	double coord, matchavg = 0;

	for (poit=0; poit<npoints; poit++)
	{
		index = nei.List[poit];
		coord = Po[index].yhex;

		if ( ( coord < reg4) && ( coord >= reg3 ) ) 
			{	
				nmatch +=1;
				matchavg += (double)Po[index].nedge/(double)Po[index].nperi;
			}

		if ( ( coord <= reg2) && ( coord > reg1 ) )
                        {
                                nmatch +=1;
                                matchavg += (double)Po[index].nedge/(double)Po[index].nperi;
                        }




	}

	if (nmatch != 0)
		{matchavg = matchavg/((double)nmatch); }
	else
		{std::cout << "There were " << nmatch << " points in the region so we move," << std::endl;
		matchavg = 1;}

	return matchavg;


}




////////////////////////////////////////////////////////////
double regionbimatchpointx(Point Po[], Nlist& nei,double reg1, double reg2)
{

	int poit=0, nmatch = 0, npoints = nei.n, index;
	double coord, matchavg = 0;

	for (poit=0; poit<npoints; poit++)
	{
		index = nei.List[poit];
		coord = Po[index].xhex;

		if ( ( coord < reg2) && ( coord >= reg1 ) ) 
			{	
				nmatch +=1;
				matchavg += Po[index].nedge/Po[index].nperi;
			}

	}

	if (nmatch != 0)
		{matchavg = matchavg/((double)nmatch); }
	else
		{std::cout << "There were " << nmatch << " cells in the region so we move," << std::endl;
		matchavg = 1;}

	return matchavg;


}


//////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
//bookmark

void CPedgealghexK(Point poit[], Point poit2[], Grid& gee,Cellprop CC[],Cellprop LL[], Hamil& engy, double Temp, double L, Nlist& nl)  
{								

//clock_t t1;

//t1 = clock();

//We get the number of points in the system

int nps = gee.NP, nall = gee.Nnum, nget = nl.n;


//std::cout << "begin function" << std::endl;

if (nget == 0 ) {std::cout << "Got 0 elements, error " << std::endl; return;}

std::uniform_int_distribution<uint32_t> uint_distAll(0,nps-1);

std::uniform_int_distribution<uint32_t> uint_distperi(0,nget-1);

std::uniform_int_distribution<uint32_t> uint_prob(1,100);

std::uniform_int_distribution<uint32_t> uint_distneihex(0,5);

engy.Hclear();		//We clear the energy


//We randomly pick a point on the grid

int ranum0 = uint_distperi(rng);

int ranum = nl.List[ranum0];

//std::cout << ranum << std::endl;

//We get the cell the point belongs to

int idc1 = poit[ranum].cnum;

//We copy the cell properties

int i=0;
int j=0;


/////////////
///
//We make a vector to get the neighbors
i = 0;
j = -1;
//We pick one of the neighbors that are of a different cell

int ic, ihe;

vector<int> neiind;

for (i=0; i<6; i++) { ic = poit[ranum].neiind[i]; ihe = poit[ic].cnum; // std::cout << "ranum neighbor " << ic << ", with cell type " << ihe << std::endl;
	if (ihe != idc1) { neiind.push_back(ic); j +=1; } 
		}

std::uniform_int_distribution<uint32_t> uint_distnei(0,j);


int ranum2 =uint_distnei(rng);

//ranum3 = neigh[ranum2]; 

int ranum3 = neiind[ranum2];

//std::cout << ranum << " , " << ranum3 << std::endl;

//We get the cell this belongs to

int idc2 = poit[ranum3].cnum; 

//And we get the original neighbor list

i = 0;
//for (i=0; i<6; i++) { ic = neigh2[i]; ihe = poit[ic].cnum; std::cout << "copy site neighbor " << ic << ", with cell type " << ihe << std::endl; }



//As well as the neighbor list for the site on the other layer
//


///////////////////////////////

//If the points belong to different cells we try to copy		
		

if  ( idc1 != idc2  )
{

int idc3 = poit2[ranum3].cnum;

i = 0;
j = 0;

//We subtract off some activity

for (j = 0; j<nps; j++) {
if (poit[j].Activ >0 ) {poit[j].Activ -= 1; } }



//The area of the cell at the change attempt changes by 1
//

CC[idc2].Ai2 = CC[idc2].Ai-1;
CC[idc1].Ai2 = CC[idc1].Ai+1;

CC[idc2].Aihex2 = CC[idc2].Aihex - sqrt(3)*L*L/2;
CC[idc1].Aihex2 = CC[idc1].Aihex + sqrt(3)*L*L/2;

//Then we introduce activity

//First is the cell extending

double gma = GMactivehex(poit,ranum);

double Hact1 = CC[idc1].lamact*gma/CC[idc1].maxact;

engy.Hactold += Hact1;
	
//Then the cell receding

gma = GMactivehex(poit,ranum3);

double Hact2 = CC[idc2].lamact*gma/CC[idc2].maxact;

engy.Hactnew += Hact2;

//std::cout << "got activity " << std::endl;

/////////////////////////////

//Now we get the change in the perimeter/adhesion/bilayer energies

int indx;

int cnuml1, cnuml2;


int Pold[nall] = {0}, Pnew[nall] = {0};
int emold1[nall] = {0}, emnew1[nall] = {0}; 
int emold2[nall] = {0}, emnew2[nall] = {0};


for (j =0; j<6; j++)
{
	indx = poit[ranum3].neiind[j];

//First the old parameters

	cnuml1 = poit[indx].cnum;

	cnuml2 = poit2[indx].cnum;

	if (cnuml1 != idc2 )

	{
		engy.Hadhold += (CC[idc2].Ji + CC[cnuml1].Ji)/2;
		Pold[idc2] += 1;
		Pold[cnuml1] += 1;

		if (cnuml2 != idc3) 
			{
				engy.Hlayold+=(CC[idc2].Jlay + LL[idc3].Jlay)/2;
				emold1[idc2] += 1;
				emold1[cnuml1] += 1;
				emold2[idc3] += 1;
				emold2[cnuml2] += 1;
			}		
		
	}


//std::cout << "got old parameters for " << j << " , " << indx << std::endl;

//Then the new parameters, cell type idc1 is now at the square

        if (cnuml1 != idc1 )

        {
                engy.Hadhnew += (CC[idc1].Ji + CC[cnuml1].Ji)/2;
                Pnew[idc1] += 1;
		Pnew[cnuml1] += 1;

                if (cnuml2 != idc3)
                        {
                                engy.Hlaynew+=(CC[idc1].Jlay + LL[idc3].Jlay)/2;
                                emnew1[idc1] += 1;
                                emnew1[cnuml1] += 1;
                                emnew2[idc3] += 1;
                                emnew2[cnuml2] += 1;

                        }

        }

//std::cout << "got new parameters for " << j << " , " << indx << std::endl;


} //End neighbor loop


//Then the area and perimeter energies. 
//
//The area of the cell at the change attempt changes by 1
//


i = 0;


for(i=0; i<nall; i++) 
{

	//The perimeter changes by the amount we tracked
//

        CC[i].Pi2 = CC[i].Pi + Pnew[i] - Pold[i];

        CC[i].edgemC2 = CC[i].edgemC + emnew1[i] - emold1[i];

	CC[i].Pihex2 = CC[i].Pihex + (Pnew[i]-Pold[i])*L/sqrt(3);


	//engy.HAold += CC[i].Asd();
	//engy.HPold += CC[i].Psd();
	engy.HAold += CC[i].Asdhex();
        engy.HPold += CC[i].Psdhex();


//	CC[i].getshapei();


	//engy.HAnew += CC[i].Asd2();
	//engy.HPnew += CC[i].Psd2();
	engy.HAnew += CC[i].Asdhex2();      
	engy.HPnew += CC[i].Psdhex2();


}

//std::cout << "Area " << engy.HAold << ", " << engy.HAnew << std::endl;
//std::cout << "Adhesion " << engy.Hadhold << " , " << engy.Hadhnew << std::endl;

//std::cout << engy.HAold << " , " << engy.HPold << " ," << engy.HAnew << " , " << engy.HPnew << std::endl; 


////////////////

//Now we do the energy check
//


double Hold = engy.getHold2();

double Hnew = engy.getHnew2();

double Bprob = exp( -(Hnew-Hold)/Temp );

//std::cout << Bprob << std::endl;

// And now we check
//
		

if (Hnew < Hold )
{

//	std::cout << "started successful copy " << std::endl;
poit[ranum3].cnum = poit[ranum].cnum;
neighpropset(poit,poit2,CC,LL,ranum3);
poit[ranum3].Activ = CC[idc1].maxact; 
i = 0;	for (i=0; i<nall; i++) 
{
	CC[i].Cellupdate(); LL[i].edgemC = LL[i].edgemC + emnew2[i]-emold2[i];
	
	if (LL[i].Pi == 0)
                {LL[i].matchratio = 0;}
        else
                {LL[i].matchratio = (double)LL[i].edgemC/(double)LL[i].Pi; }
	
}


//updatecomCC3D(gee, CC,poit[ranum3].xhex,poit[ranum3].yhex,0,CC[idc2].Aihex,idc2);
//updatecomCC3D(gee, CC,poit[ranum3].xhex,poit[ranum3].yhex,1,CC[idc1].Aihex,idc1);


j = 0;  for (j=0; j<6; j++ ) { indx = poit[ranum3].neiind[j]; 
				neighpropset(poit,poit2,CC,LL,indx);}
 

//auto k1 = find(nl.List.begin(), nl.List.end(),ranum3);

i = 0; j =0;

int nindc4, ge;

int in = 0;

for (i=0; i<nget-in; i++)
{
	nindc4 = nl.List[i];
	if (nindc4 == ranum3) { ge = nl.List[nget-1-in]; nl.List[nget-1-in] = nindc4; nl.List[i] = ge; in+=1; i-=1; }

	for (j=0; j<6; j++) {indx = poit[ranum3].neiind[j]; if (nindc4 == indx) {ge = nl.List[nget-1-in]; nl.List[nget-1-in] = nindc4; nl.List[i] = ge; in+=1; i-=1;}}

}


for (j=0; j< in; j++) { nl.List.pop_back(); }


	if (poit[ranum3].Peri == 1) { nl.List.push_back(ranum3); }
//	std::cout << "added index " << ranum3 << std::endl; ncnt = ncnt+1;}

	for (j = 0; j<6; j++) { indx = poit[ranum3].neiind[j];
	if (poit[indx].Peri == 1) { nl.List.push_back(indx); } }
//	std::cout << "added index " << nindc3 << std::endl; ncnt = ncnt+1;} }


	nl.n = nl.List.size();

//	std::cout << "end successful copy" << std::endl

}




else if (Hnew >= Hold )

	{

	double msprob = (double)uint_prob(rng)/100;

	if ( Bprob >= msprob)

	{	

//		std::cout << "begin successful copy" << std::endl;
poit[ranum3].cnum = poit[ranum].cnum;
neighpropset(poit,poit2,CC,LL,ranum3);
poit[ranum3].Activ = CC[idc1].maxact; 
i = 0; for (i=0; i<nall; i++) 
{	
	CC[i].Cellupdate(); LL[i].edgemC = LL[i].edgemC + emnew2[i]-emold2[i];
	if (LL[i].Pi == 0)
                {LL[i].matchratio = 0;}
        else
                {LL[i].matchratio = (double)LL[i].edgemC/(double)LL[i].Pi; }

}

updatecomCC3D(gee, CC,poit[ranum3].xhex,poit[ranum3].yhex,0,CC[idc2].Aihex,idc2);
updatecomCC3D(gee, CC,poit[ranum3].xhex,poit[ranum3].yhex,1,CC[idc1].Aihex,idc1);

j = 0;  for (j=0; j<6; j++ ) { indx = poit[ranum3].neiind[j]; 
		neighpropset(poit,poit2,CC,LL,indx); }
 
//auto k1 = find(nl.List.begin(), nl.List.end(),ranum3);

i = 0; j =0;

int nindc4, ge;

int in = 0;

for (i=0; i<nget-in; i++)
{
	nindc4 = nl.List[i];
	if (nindc4 == ranum3) { ge = nl.List[nget-1-in]; nl.List[nget-1-in] = nindc4; nl.List[i] = ge; in+=1; i-=1; }

	for (j=0; j<6; j++) {indx = poit[ranum3].neiind[j]; if (nindc4 == indx) {ge = nl.List[nget-1-in]; nl.List[nget-1-in] = nindc4; nl.List[i] = ge; in+=1; i-=1;}}


}

for (j=0; j< in; j++) { nl.List.pop_back(); }


	if (poit[ranum3].Peri == 1) { nl.List.push_back(ranum3); }
//	std::cout << "added index " << ranum3 << std::endl; ncnt = ncnt+1;}

	for (j = 0; j<6; j++) { indx = poit[ranum3].neiind[j];
	if (poit[indx].Peri == 1) { nl.List.push_back(indx); } }
//	std::cout << "added index " << nindc3 << std::endl; ncnt = ncnt+1;} }


	nl.n = nl.List.size();


//std::cout << "successful copy attempt" << std::endl;


	}	//Bprob check


	}	//If Hnew > Hold 


i = 0; for (i=0; i<nall; i++) { CC[i].getshapei(); CC[i].Cellupdate2(); }

engy.dH = Hnew-Hold;



} // This ends the case of two different cell types

/*
mm = 0; int nget2 = nl.n; int mnd;

for (mm = 0; mm < nget2; mm++ ) { mnd = nl.List[mm];
std::cout << "End list index " << nl.List[mm] << " perimeter " << poit[mnd].Peri << " interfaces " << poit[mnd].nperi << std::endl;}
*/

//std::cout << "successful function " << std::endl;


}	//Ends the function


/////////////////////////////////////////////////
//////////////////////////////////////////////////
//This will run the algorithm with bilayer coupling at the point level
void CPedgealgpointbend(Point poit[], Point poit2[], Grid& gee,Cellprop CC[],Cellprop LL[], Hamil& engy, double Temp, double L, Nlist& nl)  
{								

//clock_t t1;

//t1 = clock();

//We get the number of points in the system

int nps = gee.NP, nall = gee.Nnum, nget = nl.n;


//std::cout << "begin function" << std::endl;

if (nget == 0 ) {std::cout << "Got 0 elements, error " << std::endl; return;}

std::uniform_int_distribution<uint32_t> uint_distAll(0,nps-1);

std::uniform_int_distribution<uint32_t> uint_distperi(0,nget-1);

std::uniform_int_distribution<uint32_t> uint_prob(1,100);

std::uniform_int_distribution<uint32_t> uint_distneihex(0,5);

engy.Hclear();		//We clear the energy


//We randomly pick a point on the grid

int ranum0 = uint_distperi(rng);

int ranum = nl.List[ranum0];

//std::cout << ranum << std::endl;

//We get the cell the point belongs to

int idc1 = poit[ranum].cnum;

//We copy the cell properties

int i=0;
int j=0;


/////////////
///
//We make a vector to get the neighbors
i = 0;
j = -1;
//We pick one of the neighbors that are of a different cell

int ic, ihe;

vector<int> neiind;

for (i=0; i<6; i++) { ic = poit[ranum].neiind[i]; ihe = poit[ic].cnum; // std::cout << "ranum neighbor " << ic << ", with cell type " << ihe << std::endl;
	if (ihe != idc1) { neiind.push_back(ic); j +=1; } 
		}

std::uniform_int_distribution<uint32_t> uint_distnei(0,j);


int ranum2 =uint_distnei(rng);

//ranum3 = neigh[ranum2]; 

int ranum3 = neiind[ranum2];

//std::cout << ranum << " , " << ranum3 << std::endl;

//We get the cell this belongs to

int idc2 = poit[ranum3].cnum; 

//And we get the original neighbor list

i = 0;
//for (i=0; i<6; i++) { ic = neigh2[i]; ihe = poit[ic].cnum; std::cout << "copy site neighbor " << ic << ", with cell type " << ihe << std::endl; }



//As well as the neighbor list for the site on the other layer
//


///////////////////////////////

//If the points belong to different cells we try to copy		
		

if  ( idc1 != idc2  )
{

int idc3 = poit2[ranum3].cnum;

i = 0;
j = 0;

//We subtract off some activity

for (j = 0; j<nps; j++) {
if (poit[j].Activ >0 ) {poit[j].Activ -= 1; } }



//The area of the cell at the change attempt changes by 1
//

CC[idc2].Ai2 = CC[idc2].Ai-1;
CC[idc1].Ai2 = CC[idc1].Ai+1;

CC[idc2].Aihex2 = CC[idc2].Aihex - sqrt(3)*L*L/2;
CC[idc1].Aihex2 = CC[idc1].Aihex + sqrt(3)*L*L/2;

//Then we introduce activity

//First is the cell extending

double gma = GMactivehex(poit,ranum);

double Hact1 = CC[idc1].lamact*gma/CC[idc1].maxact;

engy.Hactold += Hact1;
	
//Then the cell receding

gma = GMactivehex(poit,ranum3);

double Hact2 = CC[idc2].lamact*gma/CC[idc2].maxact;

engy.Hactnew += Hact2;

//std::cout << "got activity " << std::endl;

/////////////////////////////

//Now we get the change in the perimeter/adhesion/bilayer energies

int indx;

int cnuml1, cnuml2;


int Pold[nall] = {0}, Pnew[nall] = {0};
int emold1[nall] = {0}, emnew1[nall] = {0}; 
int emold2[nall] = {0}, emnew2[nall] = {0};


double biold = poit[ranum3].Bi, biother = poit2[ranum3].Bi;

for (j =0; j<6; j++)
{
	indx = poit[ranum3].neiind[j];

//First the old parameters

	cnuml1 = poit[indx].cnum;

	cnuml2 = poit2[indx].cnum;

	if (cnuml1 != idc2 )

	{
		engy.Hadhold += (CC[idc2].Ji + CC[cnuml1].Ji)/2;
		Pold[idc2] += 1;
		Pold[cnuml1] += 1;

		if (cnuml2 != idc3) 
			{
				engy.Hlayold+=(biold+biother)/2;
				emold1[idc2] += 1;
				emold1[cnuml1] += 1;
				emold2[idc3] += 1;
				emold2[cnuml2] += 1;
			}		
		
	}


//std::cout << "got old parameters for " << j << " , " << indx << std::endl;

//Then the new parameters, cell type idc1 is now at the square


        if (cnuml1 != idc1 )

        {
                engy.Hadhnew += (CC[idc1].Ji + CC[cnuml1].Ji)/2;
                Pnew[idc1] += 1;
		Pnew[cnuml1] += 1;

                if (cnuml2 != idc3)
                        {
                               // engy.Hlaynew+=(CC[idc1].Jlay + LL[idc3].Jlay)/2;
			       engy.Hlaynew+=(biold + biother)/2;
                                emnew1[idc1] += 1;
                                emnew1[cnuml1] += 1;
                                emnew2[idc3] += 1;
                                emnew2[cnuml2] += 1;

                        }

        }

//std::cout << "got new parameters for " << j << " , " << indx << std::endl;


} //End neighbor loop


//Then the area and perimeter energies. 
//
//The area of the cell at the change attempt changes by 1
//


i = 0;


for(i=0; i<nall; i++) 
{

	//The perimeter changes by the amount we tracked
//

        CC[i].Pi2 = CC[i].Pi + Pnew[i] - Pold[i];

        CC[i].edgemC2 = CC[i].edgemC + emnew1[i] - emold1[i];

	CC[i].Pihex2 = CC[i].Pihex + (Pnew[i]-Pold[i])*L/sqrt(3);


	//engy.HAold += CC[i].Asd();
	//engy.HPold += CC[i].Psd();
	engy.HAold += CC[i].Asdhex();
        engy.HPold += CC[i].Psdhex();


//	CC[i].getshapei();


	//engy.HAnew += CC[i].Asd2();
	//engy.HPnew += CC[i].Psd2();
	engy.HAnew += CC[i].Asdhex2();      
	engy.HPnew += CC[i].Psdhex2();


}

//std::cout << "Area " << engy.HAold << ", " << engy.HAnew << std::endl;
//std::cout << "Adhesion " << engy.Hadhold << " , " << engy.Hadhnew << std::endl;

//std::cout << engy.HAold << " , " << engy.HPold << " ," << engy.HAnew << " , " << engy.HPnew << std::endl; 


////////////////

//Now we do the energy check
//


double Hold = engy.getHold2();

double Hnew = engy.getHnew2();

double Bprob = exp( -(Hnew-Hold)/Temp );

//std::cout << Bprob << std::endl;

// And now we check
//
		

if (Hnew < Hold )
{

//	std::cout << "started successful copy " << std::endl;
poit[ranum3].cnum = poit[ranum].cnum;
neighpropset(poit,poit2,CC,LL,ranum3);
poit[ranum3].Activ = CC[idc1].maxact; 
i = 0;	
for (i=0; i<nall; i++) 
{
	CC[i].Cellupdate(); LL[i].edgemC = LL[i].edgemC + emnew2[i]-emold2[i];
	
	if (LL[i].Pi == 0)
                {LL[i].matchratio = 0;}
        else
                {LL[i].matchratio = (double)LL[i].edgemC/(double)LL[i].Pi; }
	
}


updatecomCC3D(gee, CC,poit[ranum3].xhex,poit[ranum3].yhex,0,CC[idc2].Aihex,idc2);
updatecomCC3D(gee, CC,poit[ranum3].xhex,poit[ranum3].yhex,1,CC[idc1].Aihex,idc1);

j = 0;  for (j=0; j<6; j++ ) { indx = poit[ranum3].neiind[j]; 
				neighpropset(poit,poit2,CC,LL,indx);}
 

//auto k1 = find(nl.List.begin(), nl.List.end(),ranum3);

i = 0; j =0;

int nindc4, ge;

int in = 0;

for (i=0; i<nget-in; i++)
{
	nindc4 = nl.List[i];
	if (nindc4 == ranum3) { ge = nl.List[nget-1-in]; nl.List[nget-1-in] = nindc4; nl.List[i] = ge; in+=1; i-=1; }

	for (j=0; j<6; j++) {indx = poit[ranum3].neiind[j]; if (nindc4 == indx) {ge = nl.List[nget-1-in]; nl.List[nget-1-in] = nindc4; nl.List[i] = ge; in+=1; i-=1;}}

}


for (j=0; j< in; j++) { nl.List.pop_back(); }


	if (poit[ranum3].Peri == 1) { nl.List.push_back(ranum3); }
//	std::cout << "added index " << ranum3 << std::endl; ncnt = ncnt+1;}

	for (j = 0; j<6; j++) { indx = poit[ranum3].neiind[j];
	if (poit[indx].Peri == 1) { nl.List.push_back(indx); } }
//	std::cout << "added index " << nindc3 << std::endl; ncnt = ncnt+1;} }


	nl.n = nl.List.size();

//	std::cout << "end successful copy" << std::endl

}




else if (Hnew >= Hold )

	{

	double msprob = (double)uint_prob(rng)/100;

	if ( Bprob >= msprob)

	{	

//		std::cout << "begin successful copy" << std::endl;
poit[ranum3].cnum = poit[ranum].cnum;
neighpropset(poit,poit2,CC,LL,ranum3);
poit[ranum3].Activ = CC[idc1].maxact; 
i = 0; for (i=0; i<nall; i++) 
{	
	CC[i].Cellupdate(); LL[i].edgemC = LL[i].edgemC + emnew2[i]-emold2[i];
	if (LL[i].Pi == 0)
                {LL[i].matchratio = 0;}
        else
                {LL[i].matchratio = (double)LL[i].edgemC/(double)LL[i].Pi; }

}

updatecomCC3D(gee, CC,poit[ranum3].xhex,poit[ranum3].yhex,0,CC[idc2].Aihex,idc2);
updatecomCC3D(gee, CC,poit[ranum3].xhex,poit[ranum3].yhex,1,CC[idc1].Aihex,idc1);

j = 0;  for (j=0; j<6; j++ ) { indx = poit[ranum3].neiind[j]; 
		neighpropset(poit,poit2,CC,LL,indx); }
 

//auto k1 = find(nl.List.begin(), nl.List.end(),ranum3);

i = 0; j =0;

int nindc4, ge;

int in = 0;

for (i=0; i<nget-in; i++)
{
	nindc4 = nl.List[i];
	if (nindc4 == ranum3) { ge = nl.List[nget-1-in]; nl.List[nget-1-in] = nindc4; nl.List[i] = ge; in+=1; i-=1; }

	for (j=0; j<6; j++) {indx = poit[ranum3].neiind[j]; if (nindc4 == indx) {ge = nl.List[nget-1-in]; nl.List[nget-1-in] = nindc4; nl.List[i] = ge; in+=1; i-=1;}}


}

for (j=0; j< in; j++) { nl.List.pop_back(); }


	if (poit[ranum3].Peri == 1) { nl.List.push_back(ranum3); }
//	std::cout << "added index " << ranum3 << std::endl; ncnt = ncnt+1;}

	for (j = 0; j<6; j++) { indx = poit[ranum3].neiind[j];
	if (poit[indx].Peri == 1) { nl.List.push_back(indx); } }
//	std::cout << "added index " << nindc3 << std::endl; ncnt = ncnt+1;} }


	nl.n = nl.List.size();


//std::cout << "successful copy attempt" << std::endl;


	}	//Bprob check


	}	//If Hnew > Hold 


i = 0; for (i=0; i<nall; i++) { CC[i].getshapei(); CC[i].Cellupdate2(); }

engy.dH = Hnew-Hold;



} // This ends the case of two different cell types

/*
mm = 0; int nget2 = nl.n; int mnd;

for (mm = 0; mm < nget2; mm++ ) { mnd = nl.List[mm];
std::cout << "End list index " << nl.List[mm] << " perimeter " << poit[mnd].Peri << " interfaces " << poit[mnd].nperi << std::endl;}
*/

//std::cout << "successful function " << std::endl;


}	//Ends the function


///////////////////////////////////////////////////////////////
int tsfromdat(std::string& file)
{

	int nlines = 0;
	std::string words;
	std::string line;
	std::ifstream fi(file.c_str());

	while(getline(fi,line))
	{
		nlines+=1;
	}

	fi.clear();
	fi.close();

	nlines-=1;

	return nlines;

}







///////////////////////////

#endif







