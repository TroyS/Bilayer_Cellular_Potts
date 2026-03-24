//Here we setup a cell configuration using input data from a file

#include <iostream>
#include <cstring>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <memory>
#include <random>
#include <chrono> 
#include <vector>
#include "CPhexfunctions5.h"
#include <iomanip>

using std::vector;      using std::cout;        using std::cin;
using std::string;      using std::endl;        using std::ifstream;
using std::ofstream;    using std::to_string;	using std::stringstream;
using std::fixed;	using std::setprecision;


//bookmark

//begin


int main(int argc, char *argv[])


{

clock_t t1; 
t1 = clock();


initialize();


string oindx;
string infile;

if (argc > 1) {infile = argv[1];} //Input file name at the command line
if (argc > 2) { oindx = argv[2];} //Output file name at the command line

if (argc < 2) {cout << "Enter the input file name, " << endl; cin >> infile;}

if (!ifstream(infile.c_str()))
{
        cout << infile << endl;
        cout << "Could not find file." << endl;
        return 1;
}


if (argc < 3) { cout << "Enter the output file name, " << endl; cin >> oindx;}


int i=0, j=0,x=0,y=0,zlay,ntypes,nall,ncells,confchk,layeq, latchk;

ifstream inindx(infile.c_str());

string words;

inindx >> words;	inindx >> latchk;
double rad3 = sqrt(3), Latt = 1;

if (latchk != 1) { Latt = sqrt(2/rad3); } //This is the hex lattice spacing

double bend = Latt/2;           //For out of plane bending

inindx >> words; 	inindx >> confchk;

inindx >> words;	inindx >> ncells;

nall = ncells+1;

int ity = 1, ncc, nsub, nctotal=0;

vector<int> ncelv;

inindx >> words; inindx >> nsub; nctotal = nctotal+nsub; 

//Now we grab the number of grid points

inindx >> words;

for (ity=1; ity<nall; ity++)

{ inindx >> ncc; ncelv.push_back(ncc); nctotal = nctotal+ncc; }


ity = 1;

//Now we set up the region


double AAprox = nctotal*sqrt(3)*Latt*Latt/2;

int nx = ceil( sqrt(AAprox) ); 

if (nx %2 == 1 ) { nx = nx+1; }

int ny = nx;

int ns = nx*ny;

Grid gd1, gd2;

gd1.gdset2(nx, ns,nall, Latt);

gd2.gdset2(nx, ns,nall, Latt);


Cellprop celsdat1[nall], celsdat2[nall];

i = 1; 

celsdat1[0].Ainit = nsub; celsdat2[0].Ainit = nsub;

for (i=1; i<nall; i++) {celsdat1[i].Ainit = ncelv[i-1]; celsdat2[i].Ainit = ncelv[i-1];}

double eps[2];		//Adhesion between cells on the same layer
double epslay[2];	//Adhesion between cells on different layers

i = 0;

double lamsa[2];	//Lambda for area
double lamsp[2];	//Lambda for perimeter

double laac[2];		//Lambda for activity
double mac[2];		//Time length for activity

inindx >> words; inindx >> eps[0]; inindx >> eps[1]; //In plane adhesion

inindx >> words; inindx >> epslay[0]; inindx >> epslay[1]; //Bilayer adhesion

inindx >> words; inindx >> lamsa[0]; inindx >> lamsa[1]; // Area scaling

inindx >> words; inindx >> lamsp[0]; inindx >> lamsp[1]; //Perimeter scaling

inindx >> words; inindx >> laac[0]; inindx >> laac[1]; //Activity scaling

inindx >> words; inindx >> mac[0]; inindx >> mac[1]; //Activity time

if (mac[0] <= 0 ) {mac[0] == 1;} if (mac[1] <= 0 ) {mac[1] == 1;}


double Psam;		//For Perimeter

//bookmark
 
Point pnts1[ns],pnts2[ns];

int rnum, counttot, ty=1, ty2 = 1, xr, yr, hrw, movreg[ns], hx = 0;

double xcoor, ycoor;

vector<double> hex;

uniform_int_distribution<uint32_t> uint_distty(0,nall);	//Should generate random numbers between 1 and 1+ number of types
///////////////////////////

//We setup the hex lattice properties

for (hx=0; hx < ns; hx++)
{ gethxrwall(pnts1, gd1); getcoord(pnts1, gd1, hx); getcoordhex(pnts1, gd1, hx,Latt);
  gethxrwall(pnts2, gd2); getcoord(pnts2, gd2, hx); getcoordhex(pnts2, gd2, hx,Latt); 

}

//Then we setup the layers

if (nsub == 0) {
setgridVoro(pnts1, gd1, 1, 1);
setgridVoro(pnts2, gd2, 1, 2); }

if (nsub != 0) {
setgridVoro(pnts1, gd1, 0, 1);
setgridVoro(pnts2, gd2, 0, 2); }


//////////////////////////////////////////////


Cellprop celchk1[nall], celchk2[nall];

for (j=0; j<nall; j++) {celchk1[j].CPclear(); celchk2[j].CPclear(); }

int cn1, cn2;


////////////////////////

//We print the initial data to file
//

string pots =  oindx + "_init.dat";

ofstream ofile(pots.c_str());

string inithead = "npoints nx Lattchk ncells confluent \n";

string st1 = to_string(ns); string st2 = to_string(ny);

string st3 = to_string(ncells); string st4 = to_string(confchk);

inithead += st1 + " " + st2 + " " + to_string(latchk) + " " + st3 + " " + st4 + "\n";


inithead += "id x y cell_num1 cell_num2\n";


///////////////
string ot0 =   oindx + "_0.dat";

ofstream op(ot0.c_str());

string dathead = "ids type x y z vx vy vz Perimeter Bimatch Activity\n";

string regdat = "";
string datdat = "";
string datdat1 = "";
string datdat2 = "";

string st5, st6;
string dt1, dt2, dt3, dt4, dt5, dt6, dt7, dt8, dt9, dt10;

double ycor;

j = 0;

for (j = 0; j < ns; j++ )

{
	
st1 = to_string(j+1); st2 = to_string(pnts1[j].xhex);
st3 = to_string(pnts1[j].yhex); st4 = to_string(pnts1[j].cnum); 
st5 = to_string(pnts2[j].cnum);

regdat += st1 + " " + st2 + " " + st3 + " " + st4 + " " + st5 + "\n";

/////////////////

dt1 = to_string(pnts1[j].cnum); dt2 = to_string(pnts1[j].xhex);
dt3 = to_string(pnts1[j].yhex); dt4 = to_string(pnts1[j].Peri);
dt5 = to_string(pnts1[j].nedge); dt6 = to_string(pnts1[j].Activ);

datdat1 += st1 + " " + dt1 + " " + dt2 + " " + dt3 + " 0 0 0 0 " + dt4 + " " + dt5 + " " + dt6 + "\n"; 

/////////////

xcoor = pnts2[j].xhex + (nx+0.5)*Latt+2*Latt;
st6 = to_string(j+ns+1);


dt1 = to_string(pnts2[j].cnum); dt2 = to_string(xcoor);
dt3 = to_string(pnts2[j].yhex); dt4 = to_string(pnts2[j].Peri);
dt5 = to_string(pnts2[j].nedge); dt6 = to_string(pnts2[j].Activ);

datdat2 += st6 + " " + dt1 + " " + dt2 + " " + dt3 + " 0 0 0 0 " + dt4 + " " + dt5 + " " + dt6 + "\n";


}


ofile << inithead << regdat;


op << dathead << datdat1 << datdat2;

regdat.clear();
datdat.clear();
datdat1.clear();
datdat2.clear();

ofile.clear();
ofile.close();

op.clear();
op.close();


//////////


//bookmark
///////////////////////////////////////////////////////

Hamil Hlay1, Hlay2;

Hlay1.Hclear(); Hlay2.Hclear(); 

uniform_int_distribution<uint32_t> uint_distAll(0,ns-1);     //Should generate random numbers between 0 and ns-1

int e1 = 0, k = 0, rn, remchk, re, kin=0, Jco, Jinl, Joul;

double HH;

vector<int> nlist1;
vector<int> nlist2;

Hlay1.Hclear(); Hlay2.Hclear();

//We count the number of grid points in each cell
ity = 0;

int ts = 0, tn=0, tn2=0, tsteps;

//We should let the temperature be defined by the user

double T;

inindx >> words; inindx  >> T;

inindx >> words; inindx >> tsteps; 

string tprint, ot, lays;

int ggs = 0, ccount1, ccount2, ui = 1;

vector<double> coms;

string opro = "props"; 

//bookmark
///////////////////////////////////////////////////////////////
//We get the area

int bb = 0;

for (bb=0; bb < nall; bb++ )
{ celsdat1[bb].kapp = bend; celsdat2[bb].kapp = bend;}

bb = 0;

Nlist nei1, nei2;               //This is for the neighbor list

//Then we run though both layers to get the original cell areas and perimeters, and to set up the neighbor list

setncprops(pnts1,pnts2,celsdat1,celsdat2,gd1,gd2,nei1,nei2,Latt);

getcomCC3D(pnts1,gd1,celsdat1,1);
getcomCC3D(pnts2,gd2,celsdat2,1);
hexorder(pnts1,celsdat1,gd1,1);
hexorder(pnts2,celsdat2,gd2,1);

int l1 = 0, l2 = 0;

i = 0;

for (i=0; i<nall; i++ ) { 
 cout << "The hex area of cell " << i << " in layer 1 is " << celsdat1[i].Aihex << " , the hex perimeter is " << celsdat1[i].Pihex << endl;

cout << endl;

cout << "The hex area of cell " << i << " in layer 2 is " << celsdat2[i].Aihex << ", the hex perimeter is " << celsdat2[i].Pihex << endl;

cout << endl;
	

}


cout << "Now we start the initialization, " << endl;
cout << endl;

double Perii1[nall] = {0}, Perii2[nall] = {0};

double areaf1[nall], areaf2[nall], Aa;

inindx >> words;

for (i=0; i<ncells; i++)

{ inindx >> areaf1[i]; }

inindx >> words;

for (i=0; i<ncells; i++)

{ inindx >> areaf2[i]; }


/////////////

inindx >> words;

for (i=0; i<ncells; i++)

{inindx >> Perii1[i]; }

inindx >> words;

for (i=0; i<ncells; i++)

{inindx >> Perii2[i]; }


/////////////////////////////

//We use cell 0 for the substrate 
//
//
 
celsdat1[0].CPsetall(0, 0, 0, 0, eps[0], 0, laac[0], mac[0]);
celsdat1[0].lay = 1;



celsdat2[0].CPsetall(0, 0, 0, 0, eps[0], 0, laac[0], mac[0]);
celsdat2[0].lay = 2;

celsdat1[0].A0 = celsdat1[0].Ai; celsdat2[0].A0 = celsdat2[0].Ai;


i = 1;

for (i=1; i<nall; i++ )
{
	
celsdat1[i].lay = 1; 
celsdat1[i].CPsetall(lamsa[0], areaf1[i-1], lamsp[0], Perii1[i-1], eps[0], epslay[0], laac[0], mac[0]);

celsdat2[i].lay = 2; 
celsdat2[i].CPsetall(lamsa[1], areaf2[i-1], lamsp[1], Perii2[i-1], eps[1], epslay[1], laac[1], mac[1]); 
 

}


/////////////////
//Let's make a file for the matching edges vs timesteps
//
string edmtpr =  oindx + "edgematch.dat";

ofstream epr(edmtpr.c_str());

string ematch = "J lambdaA A0lay1 A0lay2 lambdaP P0lay1 P0lay2 lambdaBi\n";


string et1 = to_string(eps[0]); string et2 = to_string(lamsa[0]); 
string et3 = to_string(areaf1[0]); string et4 = to_string(areaf2[1]);
string et5 = to_string(lamsp[0]); string et6 = to_string(Perii1[1]);
string et7 = to_string(Perii2[1]); string et8 = to_string(epslay[0]);

ematch += et1 + " " + et2 + " " + et3 + " " + et4 + " " + et5 + " " + et6 + " " + et7 + " " + et8 + "\n";


ematch += "Timestep Totperim Totedgematch Edgematchratio Avghexorder1 Avghexorder2 dElay1 dElay2\n";


double Perimavg = 0, Edgematavg = 0, Matchperiratio = 0, Horder1, Horder2;
int freq;
inindx >> words >> freq;

dathead = "ids type x y z vx vy vz Perimeter Bimatch Activity\n";

/////////////////////////////////
//bookmark
//Here we start the CPM

for (ts = 0; ts < tsteps; ts++)

{						// The MMC step

Perimavg = 0;	Edgematavg = 0;	Matchperiratio = 0;


	for (tn = 0; tn < ns; tn++ )


	{ 				// The pixel loop

CPedgealghexK(pnts1, pnts2, gd1, gd2, celsdat1,celsdat2,  Hlay1, Hlay2, T,Latt, nei1);


CPedgealghexK(pnts2, pnts1, gd2, gd1, celsdat2,celsdat1,  Hlay2, Hlay1, T,Latt, nei2);

	}		//End loop for algorithm

//We print for paraview
//
	
	if  ( ts % freq == 0 ) 

	{

		tprint = to_string(ts+freq);

		ot =  oindx + "_" +tprint + ".dat";

		ggs = 0;

		ofstream ogs(ot.c_str());

		datdat1 = "";
		datdat2 = "";
	

		for (ggs=0; ggs < ns; ggs++ )

		{


			dt1 = to_string(ggs+1); dt2 = to_string(pnts1[ggs].cnum);
			dt3 = to_string(pnts1[ggs].xhex); dt4 =to_string(pnts1[ggs].yhex);
			dt5 = to_string(pnts1[ggs].zhex); dt6=to_string(pnts1[ggs].Peri);
			dt7 = to_string(pnts1[ggs].nedge);dt8=to_string(pnts1[ggs].Activ);

			datdat1+= dt1 + " " + dt2 + " " + dt3 + " " + dt4 + " " + dt5 			+" " + " 0 0 0 " + dt6 + " " + dt7 + " " + dt8 + "\n";



////////////////////////////

			xcoor = pnts2[ggs].xhex + (nx+0.5)*Latt + 2*Latt;

			dt1 = to_string(ggs+1+ns); dt2 = to_string(pnts2[ggs].cnum);
                	dt3 = to_string(xcoor); dt4 =to_string(pnts2[ggs].yhex);
                	dt5 = to_string(pnts2[ggs].zhex); dt6=to_string(pnts2[ggs].Peri);
                	dt7 = to_string(pnts2[ggs].nedge);dt8=to_string(pnts2[ggs].Activ);

			datdat2+= dt1 + " " + dt2 + " " + dt3 + " " + dt4 + " " + dt5 			+ " " + " 0 0 0 " + dt6 + " " + dt7 + " " + dt8 + "\n";


        	}	//End loop over sites for printing

	ogs << dathead << datdat1 << datdat2;

        ogs.clear();
        ogs.close();

////////////////


	}		//Paraview loop

		
getcomCC3D(pnts1,gd1,celsdat1,1);
getcomCC3D(pnts2,gd2,celsdat2,1);
hexorder(pnts1,celsdat1,gd1,1);
hexorder(pnts2,celsdat2,gd2,1);
		
Horder1 = meanhexorder(celsdat1, gd1);
Horder2 = meanhexorder(celsdat2, gd2);
		
ui = 0; 
		
for(ui=0; ui<nall; ui++)
{
	Perimavg = Perimavg + celsdat1[ui].Pi + celsdat2[ui].Pi;
        Edgematavg = Edgematavg + celsdat1[ui].edgemC + celsdat2[ui].edgemC;
}

Perimavg = Perimavg/(2);
Edgematavg = Edgematavg/(2);

Matchperiratio = Edgematavg/Perimavg;


et1 = to_string(ts); et2 = to_string(Perimavg);
et3 = to_string(Edgematavg); et4 = to_string(Matchperiratio);
et5 = to_string(Horder1); et6 = to_string(Horder2);
et7 = to_string(Hlay1.dH); et8 = to_string(Hlay2.dH);

ematch+= et1 + " " + et2 + " " + et3 + " " + et4 + " " + et5 + " " + et6 + " " + et7 + " " + et8 + "\n"; 


}


epr << ematch;

epr.clear();
epr.close();

ematch.clear();

//////////////////////////////////////////


//We print to file
//

string potsout =  oindx + "region.dat";

ofstream ofile2(potsout.c_str());

regdat = "";

j = 0;

ccount1 = 0;

ccount2 = 0;

for (j = 0; j < ns; j++ )

	{
	
	regdat += to_string(j) + " " + to_string(pnts1[j].xhex) + " " + to_string(pnts1[j].yhex) + " " + to_string(pnts1[j].cnum) + " " + to_string(pnts2[j].cnum)+ "\n";

	}


ofile2 << inithead << regdat;

regdat.clear();

ofile2.clear();
ofile2.close();




//////////////////////////////////////////////////////////


//We clear the cell data
//bookmark


vector<int> check, check2;

ity = 0;

int medge1 = 0, medge2 = 0;

for (ity=0; ity<ns; ity++ )
{medge1 = medge1 + pnts1[ity].nedge; medge2 = medge2 + pnts2[ity].nedge; }


ity = 0;

int ity2 = 0;


//Now we allow changes for the desired area

i = 0;

string ostat =  oindx + "stats.txt";

ofstream sta(ostat.c_str());

for (i=0; i<nall; i++ ) 
{ 
cout << "The hex area of cell " << i << " in layer 1 is " << celsdat1[i].Aihex << " , the hex perimeter is " << celsdat1[i].Pihex <<  
" , and the hexatic order is " << celsdat1[i].hexaorder << endl;	

sta << "Cell" << i << "lay1" << " " << "Area" << " " << celsdat1[i].Aihex << " " << "Perimeter" << " " << celsdat1[i].Pihex << " " << "hexorder"  << " " << celsdat1[i].hexaorder << " " << "xcom" << " " << celsdat1[i].xcom << " " << "ycom" << " " << celsdat1[i].ycom <<endl;

cout << endl;
//sta << endl;

cout << "The hex area of cell " << i << " in layer 2 is " << celsdat2[i].Aihex << ", the hex perimeter is " << celsdat2[i].Pihex << " , and the hexatic order is " << celsdat2[i].hexaorder << endl;
	
sta << "Cell" << i << "lay2" << " " << "Area" << " " << celsdat2[i].Aihex << " " << "Perimeter" << " " << celsdat2[i].Pihex << " " << "hexorder"  << " " << celsdat2[i].hexaorder << " " << "xcom" << " " << celsdat2[i].xcom << " " << "ycom" << " " << celsdat2[i].ycom <<endl;



//sta << endl;

cout << endl;


}

///////////////


cout << "The average number of matching edges is " << (medge1+medge2)/(2*(double)nall) << endl;


sta.clear();
sta.close();

i = 1;


t1 = clock()-t1;

cout << t1/(60*CLOCKS_PER_SEC) << " minutes " << t1/(CLOCKS_PER_SEC) << " seconds. " << endl;




}
