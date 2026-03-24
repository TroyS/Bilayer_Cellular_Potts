//Here we run a parameter sweep on bilayer coupling strength

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
#include <cmath>

using std::vector;      using std::cout;        using std::cin;
using std::string;      using std::endl;        using std::ifstream;
using std::ofstream;    using std::to_string;	using std::stringstream;
using std::fixed;	using std::setprecision;


//bookmark

//begin


int main(int argc, char *argv[])


{

clock_t t1; 
clock_t t2;
t1 = clock();


initialize();

// We get the input file

string oindx;
string pfile;
string infile;

if (argc > 3) {infile = argv[1]; pfile = argv[2]; oindx = argv[3];} 
// For the system setup and parameters


if (argc < 4) 
{
	cout << "Enter the input file name, " << endl;
	cin >> infile; cout << "Enter the parameter file name, " << endl; 
	cin >> pfile; cout << "Enter the output file name, " << endl;
	cin >> oindx;
}

if (!ifstream(infile.c_str()))
{
        cout << infile << endl;
        cout << "Could not find file." << endl;
        return 1;
}


if (!ifstream(pfile.c_str()))
{
        cout << pfile << endl;
        cout << "Could not find parameter file." << endl;
        return 1;
}


if (argc < 3) { cout << "Enter the output file name, " << endl;
cin >>  oindx;}


//////////////////////////////////


int i = 0,j = 0, nall,ncells,nx,ny,ns,la,hw,latchk,confchk;

string words;

double Da,Ya,Latt = 1,rad3 = sqrt(3);


ifstream fi(infile.c_str());

fi >> words >> words >> words >> words >> words; //Get the header
fi >> ns >> nx >> latchk >> ncells >> confchk;
fi >> words >> words >> words >> words >> words; //Get the second header

if (latchk != 1) { Latt = sqrt(2/rad3); } //This is the hex lattice spacing

nall = ncells+1;

double bend = Latt/2;

ny = ns/nx;

Grid gd1, gd2;

gd1.gdset3(nx, ny,nall,Latt);

gd2.gdset3(nx, ny,nall,Latt);


Cellprop celstart1[nall], celstart2[nall];

Point pstart1[ns],pstart2[ns];

int xcord, ycord;


for (i=0; i< ns; i++)
{

        xcord = (i%nx)+1; ycord = floor(i/nx) + 1;


        fi >> la; pstart1[i].indx = la; pstart2[i].indx = la;

        pstart1[i].xcart = xcord; pstart2[i].xcart = xcord;

        pstart1[i].ycart = ycord; pstart2[i].ycart = ycord;

        fi >> Da; pstart1[i].xhex = Da; pstart2[i].xhex = Da;

        fi >> Ya; pstart1[i].yhex = Ya; pstart2[i].yhex = Ya;

         hw = ycord % 2;

        if (hw==1) {pstart1[i].hexrw = 1; pstart2[i].hexrw = 1; }

        if (hw==0) {pstart1[i].hexrw = 0; pstart2[i].hexrw = 0; }


        fi >> pstart1[i].cnum;    fi >> pstart2[i].cnum;

        pstart1[i].zhex = 0; pstart2[i].zhex = 0;

}


////////////////

//Then we get the props

ifstream inindx(pfile.c_str());

int skip;

inindx >> words >> skip;
inindx >> words >> skip;
inindx >> words >> skip;
inindx >> words >> skip;
inindx >> words;

for (i=0; i<ncells; i++) {inindx >> skip; }

i = 1; 

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

if (mac[0] <= 0 ) {mac[0] = 1;} if (mac[1] <= 0 ) {mac[1] = 1;}

int ity = 0, ts = 0, tn = 0, tsteps;

double T, Peris1[nall] = {0}, Peris2[nall] = {0}, areas1[nall], areas2[nall];

inindx >> words; inindx  >> T;

inindx >> words; inindx >> tsteps; 


inindx >> words;
 
for (i=0; i<ncells; i++)

{ inindx >> areas1[i]; }

inindx >> words;

for (i=0; i<ncells; i++)

{ inindx >> areas2[i]; }


/////////////

inindx >> words;

for (i=0; i<ncells; i++)

{inindx >> Peris1[i]; }

inindx >> words;

for (i=0; i<ncells; i++)

{inindx >> Peris2[i]; }

int freq;

inindx >> words >> freq;


//Then we get the parameters for the sweep

double Biend; inindx >> words >> Biend; inindx >> Biend;

double Bifreq; inindx >> words >> Bifreq;

double xcoor;

Hamil Hlay1, Hlay2;
vector<int> nlist1,nlist2;

string tprint,ot, lays;

int ggs = 0, ui =1, bb=0;

string opro = "props";


for (bb=0; bb < nall; bb++ )
{ celstart1[bb].kapp = bend; celstart2[bb].kapp = bend;}

Nlist nstart1, nstart2, nei1, nei2;              //This is for the neighbor list

//Then we run though both layers to get the original cell areas and perimeters, and to set up the neighbor list

setncprops(pstart1,pstart2,celstart1,celstart2,gd1,gd2,nstart1,nstart2,Latt);

getcomCC3D(pstart1,gd1,celstart1);
getcomCC3D(pstart2,gd2,celstart2);
hexorder(pstart1,celstart1,gd1);
hexorder(pstart2,celstart2,gd2);


i = 0;

for (i=0; i<nall; i++ ) { 
 cout << "The hex area of cell " << i << " in layer 1 is " << celstart1[i].Aihex << " , the hex perimeter is " << celstart1[i].Pihex <<  " , and the hexatic order is " << celstart1[i].hexaorder << endl;
	 
cout << endl;

cout << "The hex area of cell " << i << " in layer 2 is " << celstart2[i].Aihex << ", the hex perimeter is " << celstart2[i].Pihex << " , and the hexatic order is " << celstart2[i].hexaorder << endl;

cout << endl;


}


////////////////////////////////////////////////
//We'll print the inital setup to file
string ot0;

string dathead = "ids type x y z vx vy vz Perimeter Bimatch Activity\n";

string datinit1 = "";
string datinit2 = "";
string datdat1 = "";
string datdat2 = "";

string dt1, dt2, dt3, dt4, dt5, dt6, dt7, dt8, dt9, dt10;

string st1, st6, ostat;

j = 0;

for (j = 0; j < ns; j++ )

{

st1 = to_string(j+1);

dt1 = to_string(pstart1[j].cnum); dt2 = to_string(pstart1[j].xhex);
dt3 = to_string(pstart1[j].yhex); dt4 = to_string(pstart1[j].Peri);
dt5 = to_string(pstart1[j].nedge); dt6 = to_string(pstart1[j].Activ);

datinit1 += st1 + " " + dt1 + " " + dt2 + " " + dt3 + " 0 0 0 0 " + dt4 + " " + dt5 + " " + dt6 + "\n";

//////////////////

xcoor = pstart2[j].xhex + (nx+0.5)*Latt+2*Latt;

dt7 = to_string(j+ns+1);


dt1 = to_string(pstart2[j].cnum); dt2 = to_string(xcoor);
dt3 = to_string(pstart2[j].yhex); dt4 = to_string(pstart2[j].Peri);
dt5 = to_string(pstart2[j].nedge); dt6 = to_string(pstart2[j].Activ);

datinit2 += dt7 + " " + dt1 + " " + dt2 + " " + dt3 + " 0 0 0 0 " + dt4 + " " + dt5 + " " + dt6 + "\n";

	}

vector<int> check, check2;

int medge1 = 0, medge2 = 0;

double Perii1[nall] = {0}, Perii2[nall] = {0},areaf1[nall], areaf2[nall];
 
for (i=0; i<ncells; i++)

{ areaf1[i] = areas1[i]; }

for (i=0; i<ncells; i++)

{ areaf2[i] = areas2[i]; }


/////////////

for (i=0; i<ncells; i++)

{ Perii1[i] = Peris1[i]; }

for (i=0; i<ncells; i++)

{ Perii2[i] = Peris2[i]; }





//bookmark
////////////////////////////////////////
////////////////////////////////////////////////bookmark
//Now, we start the bilayer parameter sweep

int swp = 0;

//double dswp = (Biend - epslay[0])/Bifreq;

int nswp = abs(ceil((Biend - epslay[0])/Bifreq));

double Bipara,Perimavg,Edgematavg, Matchperiratio,Horder1, Horder2;

string Bilabel,para,edmtpr,ematch; 

string et1, et2, et3, et4, et5, et6, et7, et8;

string regdat;
string inithead = "npoints nx Lattchk ncells confluent \n";
string st2, st3, st4;

Point pnts1[ns], pnts2[ns];

Cellprop celsdat1[nall], celsdat2[nall];

for (swp = 0; swp < nswp; swp++)

{

nei1 = nstart1; nei2 = nstart2;

Bipara = epslay[0]+Bifreq*swp;
stringstream stream;
stream << fixed << setprecision(1) << std::abs(Bipara);
Bilabel = stream.str();
para = "Bi_" + Bilabel;

///////////////////////////////////
////////////////////////
//We'll print the inital setup to file
ot0 =  oindx + para + "_0.dat";

ofstream op(ot0.c_str());
op << dathead << datinit1 << datinit2;

op.clear();
op.close();

///////////////////////////////////


Hlay1.Hclear(); Hlay2.Hclear(); 

for (i=0; i<ns; i++) {pnts1[i] = pstart1[i]; pnts2[i] = pstart2[i]; }

for (j=0; j<nall; j++) {celsdat1[j] =celstart1[j]; celsdat2[j] = celstart2[j];}

bb = 0;

///////////////////////////////////
/////////////////////////////

//We use cell 0 for the substrate 
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
celsdat1[i].CPsetall(lamsa[0], areaf1[i-1], lamsp[0], Perii1[i-1], eps[0], Bipara, laac[0], mac[0]);

celsdat2[i].lay = 2; 
celsdat2[i].CPsetall(lamsa[1], areaf2[i-1], lamsp[1], Perii2[i-1], eps[1], Bipara, laac[1], mac[1]); 
 

}



/////////////////
//Let's make a file for the matching edges vs timesteps
//
edmtpr =  oindx + para + "edgematch.dat";

ofstream epr(edmtpr.c_str());

ematch = "J lambdaA A0lay1 A0lay2 lambdaP P0lay1 P0lay2 lambdaBi\n";

et1 = to_string(eps[0]); et2 = to_string(lamsa[0]);// et2 = to_string(eps[1]); 
et3 = to_string(areaf1[0]); et4 = to_string(areaf2[1]);
et5 = to_string(lamsp[0]); et6 = to_string(Perii1[1]);
et7 = to_string(Perii2[1]); et8 = Bilabel;  

ematch += et1 + " " + et2 + " " + et3 + " " + et4 + " " + et5 + " " + et6 + " " + et7 + " " + et8 + "\n";

ematch += "Timestep Totperim Totedgematch Edgematchratio Avghexorder1 Avghexorder2 dElay1 dElay2\n";


Perimavg = 0;
Edgematavg = 0;
Matchperiratio = 0;

datdat1 = "";
datdat2 = "";
regdat = "";
inithead = "npoints nx Lattchk ncells confluent \n";
st1 = to_string(ns); st2 = to_string(nx);

st3 = to_string(ncells); st4 = to_string(confchk);

inithead += st1 + " " + st2 + " " + to_string(latchk) + " " + st3 + " " + st4 + "\n";

cout << "Now we start the run for  " << para << " ," << endl;
cout << endl;

ts = 0; tn = 0;

for (ts = 0; ts < tsteps; ts++)

{						// The MMC step

Perimavg = 0; Edgematavg = 0; Matchperiratio = 0;


	for (tn = 0; tn < ns; tn++ )


	{ 				// The pixel loop

CPedgealghexK(pnts1, pnts2, gd1, celsdat1,celsdat2, Hlay1, T,Latt, nei1);

CPedgealghexK(pnts2, pnts1, gd2, celsdat2,celsdat1, Hlay2, T,Latt, nei2);

	}		//End loop for algorithm


//We print for paraview


	if  ( ts % freq == 0 ) 

	{


		tprint = to_string(ts+freq);

		ot =  oindx + para + "_" +tprint + ".dat";

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


getcomCC3D(pnts1,gd1,celsdat1);
getcomCC3D(pnts2,gd2,celsdat2);

hexorder(pnts1,celsdat1,gd1);
hexorder(pnts2,celsdat2,gd2);
		

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

stream.str("");


//////////////////////////////////////////


//We print to file
//

string potsout =  oindx + para + "region.dat";

ofstream ofile2(potsout.c_str());

regdat = "";
j = 0;
;

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

ity = 0;

medge1 = 0; medge2 = 0;

for (ity=0; ity<ns; ity++ )
{medge1 = medge1 + pnts1[ity].nedge; medge2 = medge2 + pnts2[ity].nedge; }


ity = 0;

i = 0;

ostat =  oindx + para + "stats.txt";

ofstream sta(ostat.c_str());

for (i=0; i<nall; i++ ) 
{
	cout << "The hex area of cell " << i << " in layer 1 is " << celsdat1[i].Aihex << " , the hex perimeter is " << celsdat1[i].Pihex << " , and the hexatic order is " << celsdat1[i].hexaorder << endl;

sta << "Cell" << i << "lay1" << " " << "Area" << " " << celsdat1[i].Aihex << " " << "Perimeter" << " " << celsdat1[i].Pihex << " " << "hexorder"  << " " << celsdat1[i].hexaorder << " " << "xcom" << celsdat1[i].xcom << " " << "ycom" << celsdat1[i].ycom <<endl;

cout << endl;
//sta << endl;

cout << "The hex area of cell " << i << " in layer 2 is " << celsdat2[i].Aihex << ", the hex perimeter is " << celsdat2[i].Pihex << " , and the hexatic order is " << celsdat2[i].hexaorder << endl;

sta << "Cell" << i << "lay2" << " " << "Area" << " " << celsdat2[i].Aihex << " " << "Perimeter" << " " << celsdat2[i].Pihex << " " << "hexorder"  << " " << celsdat2[i].hexaorder << " " << "xcom" << celsdat2[i].xcom << " " << "ycom" << celsdat2[i].ycom <<endl;

cout << endl;


}

cout << "The average number of matching edges is " << (medge1+medge2)/(2*(double)nall) << endl;


sta.clear();
sta.close();

i = 1;


t2 = clock()-t1;

cout << t2/(60*CLOCKS_PER_SEC) << " minutes " << t2/(CLOCKS_PER_SEC) << " seconds. " << endl;

} //End parameter sweep


}
