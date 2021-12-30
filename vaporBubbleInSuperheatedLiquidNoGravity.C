/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2010-08-02 Eelco van Vliet: 1st public version of vaporBubbleInSuperheatedLiquidNoGravity:
  http://www.cfd-online.com/Forums/openfoam-solving/66705-wallheatflux-bc-not-constant-after-restart.html#post269812

2012-05-21 Eelco van Vliet:
  Quoting: http://www.cfd-online.com/Forums/openfoam-post-processing/101972-wallheatflux-utility-incompressible-case.html#post362191
  «modified the standard wallHeatflux utility which comes default with OF into
  a version for incompressible flows. Also removed a bug out of the code.»

2012-06-26 Eelco van Vliet:
  Quoting: http://www.cfd-online.com/Forums/openfoam-post-processing/101972-wallheatflux-utility-incompressible-case.html#post368330
  «p is now not required anymore.»

2014-06-22: Bruno Santos: Adapted to OpenFOAM 2.2.x.

-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    vaporBubbleInSuperheatedLiquidNoGravity

Description
    Calculates and writes the the analytical solution for temperature
	in case of growing of vapor bubble in a superheated liquid under
	zero gravity conditions.

References:
@article{SATO2013127,
title = {A sharp-interface phase change model for a mass-conservative interface tracking method},
journal = {Journal of Computational Physics},
volume = {249},
pages = {127-161},
year = {2013},
issn = {0021-9991},
doi = {https://doi.org/10.1016/j.jcp.2013.04.035},
url = {https://www.sciencedirect.com/science/article/pii/S0021999113003197},
author = {Yohei Sato and Bojan Ničeno},
keywords = {Phase change, Interface tracking, Conservative method, Mass transfer, Color function},
abstract = {A new phase-change model has been developed for a mass-conservative interface tracking method. The mass transfer rate is directly calculated from the heat flux at the liquid–vapor interface, and the phase change takes place only in the cells which include this interface. As a consequence of the sharpness of the mass transfer rate distribution, the velocity jump across the interface can be captured, and high accuracy can be maintained. The method has been implemented in an incompressible Navier–Stokes equations solver employing a projection method based on a staggered finite-volume algorithm on Cartesian grids. The model has been verified for one-dimensional phase-change problems and a three-dimensional simulation of a growing vapor bubble in a superheated liquid under zero gravity condition. The computed results agree with theoretical solutions, and the accuracy of the model is confirmed to be of second-order in space using a grid refinement study. A three-dimensional simulation of a rising vapor bubble in a superheated liquid under gravity has been performed as a validation case, and good agreement with experimental data is obtained for the bubble growth rate. As a demonstration of the applicability of the method to engineering problems, a nucleate boiling simulation is presented with a comparison to experimental data. Good agreement is obtained for the bubble shapes and the bubble departure period. In all the simulation cases, strict mass conservation is satisfied.}
}
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
// modified from  wallHeatFlux
#include "singlePhaseTransportModel.H"
#include "phaseChangeTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "turbulenceModel.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Calculates analtical value of bubble radius
dimensionedScalar bubbleRadius
       (
	       const dimensionedScalar beta, 
		   const dimensionedScalar liquidThermCond, 
		   const dimensionedScalar liquidRho, 
		   const dimensionedScalar liquidCp, 
		   const scalar time
	   )
{
	return 2*beta*Foam::sqrt(liquidThermCond*time/liquidRho/liquidCp);
}

// Calculates integral (Eq. (34) or Eq. (35))
dimensionedScalar calcIntegral
(
	const label N,                        // number of rectangles
	const dimensionedScalar beta,         // beta_g
	const dimensionedScalar R,            // initial bubble radius
	const dimensionedScalar r,            // radius 
	const dimensionedScalar rho1,         // liquid density
	const dimensionedScalar rho2          // vapor density
)
{
	// Number of nodes
	const label nodes = N+1;
	
	// Thickness of a rectangular
	const scalar dksi = ( 1.0-(1-(R/r).value()) )/N;
	scalar ksi = 1-(R/r).value();
	dimensionedScalar integral("integral0", dimless, 0);

	for (int i = 0; i<nodes-1; i++)
	{
		integral += Foam::exp
		(
			-Foam::pow(beta, 2)*(Foam::pow(1-ksi, -2) - 2*(1 - rho2/rho1)*ksi - 1)
		)*dksi;
		ksi += dksi;
	}

	return integral;
}

int main(int argc, char *argv[])
{

    timeSelector::addOptions();
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"
	#include "readTransportProperties.H"
    #include "createFields.H"

	const label maxIterNr = 20000;
	label iter = 0;
	const scalar tol = 1e-12;
    
	// LHS of Eq. (34)
	const dimensionedScalar LHS = rho1*cp1*(Tinf-TSat)/rho2/(hEvap+(cp1-cp2)*(Tinf-TSat));
	
	// Integral in Eq. (34)
	dimensionedScalar RHS("RHS", dimless, 0.0);
	dimensionedScalar beta_g_prev = 0;               // guess of starting value

    Info<< endl << "Obliczam beta_g... " << endl;
	do 
	{
		beta_g_prev = beta_g;		
		RHS = calcIntegral(N, beta_g, R, R, rho1, rho2);
		beta_g = Foam::sqrt(LHS/2.0/RHS);
		iter++;
	} while ( mag(beta_g.value() - beta_g_prev.value()) > tol && iter < maxIterNr);

    Info<< endl << "Liczba iteracji: " << iter << endl;
	Info<< "beta_g = " << beta_g.value() << endl;

    Info<< endl << "Calculating analytical temperature distribution... " << endl;
	OFstream IFfileT("bubbleTemperature.txt");
	IFfileT << "Radius [m]\t" << "Analytical temperature [K]" << endl;
    dimensionedScalar Tempr0("Tempr0", dimTemperature, 0);
	List<dimensionedScalar> Tempr(Rsegments+1, Tempr0.value());
	dimensionedScalar dr = (Rend-R)/Rsegments;	
	const dimensionedScalar B = rho2*(hEvap+(cp1-cp2)*(Tinf-TSat))/rho1/cp1;

	List<dimensionedScalar> r(Rsegments+2, R.value());
	for (int rstep=0; rstep<Rsegments+1; rstep++)
	{
		dimensionedScalar integral("integral", dimless, 0);
		integral = calcIntegral(N, beta_g, R, r[rstep], rho1, rho2);
		Tempr[rstep] = Tinf.value() - 2*beta_g*beta_g*B.value()*integral;
		IFfileT << r[rstep].value() << "\t" << Tempr[rstep].value() << endl;
		r[rstep+1] = r[rstep] + dr.value();
	}
	Info<< "\nSaving the results to bubbleTemperature.txt\n" << endl;


	if (calcTiniAnalytical)
	{
		Info<< endl << "Creating Tini file in 0 directory with analytical temperature distribution... " << endl;
		// Creating initial conditions for temperature based on analytical solution
    	volScalarField radius
    	(
    	    IOobject
    	    (
				"radius",
    	        runTime.timeName(),
    	        mesh,
    	        IOobject::NO_READ,
    	        IOobject::NO_WRITE
    	    ),
    	    mesh,
			dimensionedScalar("r0", dimLength, 0)
    	);

		const volVectorField cells = mesh.C();

		forAll(cells, celli)
		{
			radius[celli] = Foam::sqrt
			(
				Foam::pow(cells[celli].component(0) - Rorigin.component(0).value(), 2)
		      + Foam::pow(cells[celli].component(1) - Rorigin.component(1).value(), 2)
		      + Foam::pow(cells[celli].component(2) - Rorigin.component(2).value(), 2)
			);
		}

		volScalarField Tini
		(
		    IOobject
		    (
		        "Tini",
		        runTime.timeName(),
		        mesh,
		        IOobject::NO_READ,
		        IOobject::AUTO_WRITE
		    ),
		    T
		);

		//// For linear initialization
		//const scalar b = (Tinf.value()-TSat.value()*Rend.value()/R.value())/(1.0-Rend.value()/R.value());
		//const scalar a = (TSat.value()-b)/R.value();
		const scalar conv = 1e-5;
		Tini = Tinf;
		forAll(Tini, celli)
		{
			if (radius[celli] < R.value()) 
			{
				Tini[celli] = TSat.value();
			}
		}
		forAll(Tini, celli)
		{
			for (int rstep=0; rstep<Rsegments+1; rstep++)
			{
				if (mag(Tempr[rstep]-Tinf.value()).value() < conv) continue; 
				if (radius[celli] >= r[rstep].value() 
				 && radius[celli] <  r[rstep+1].value() )
				{
					Tini[celli] = Tempr[rstep].value();
					//Tini[celli] = a*r[rstep].value() + b;
				}
			}
		}
		Tini.write();
	}

	Info<< "Calculating numerical and analytical radius... \n" << endl;
	OFstream IFfile("bubbleRadius.txt");
	IFfile << "Time [s]\t" << "Numerical [m]\t" << "Analytical [m]\t" << "Error [%]" << endl;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< endl << "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

		Info<< "Reading field alpha." << mixture->phase1Name() << endl;
    	volScalarField alphal
    	(
    	    IOobject
    	    (
				"alpha." + mixture->phase1Name(),
    	        runTime.timeName(),
    	        mesh,
    	        IOobject::MUST_READ,
    	        IOobject::NO_WRITE
    	    ),
    	    mesh
    	);

        gradAlphal=fvc::grad(alphal);

        Info<< endl;

        dimensionedScalar radius("radius", dimLength, 0.0);
	    radius = sum(mag(gradAlphal)*mesh.C().component(gradAlphaCalcDir))/sum(mag(gradAlphal));

		Info<<"Numerical value of bubble radius for time " 
			<< runTime.timeName() 
			<< " is equal to: " 
			<< radius.value() 
			<< " [m]"
			<< " (error: "
			<< mag(radius.value() - bubbleRadius(beta_g, k1, rho1, cp1, runTime.value()).value())/
					mag(bubbleRadius(beta_g, k1, rho1, cp1, runTime.value()).value() + VSMALL)*100
			<< "%)"
			<< endl;

	    Info<< "\nSaving the results to bubbleRadius.txt\n" << endl;

		IFfile << runTime.timeName() 
	  	     << "\t" 
			 << radius.value() 
	  	     << "\t" 
	  		 << bubbleRadius(beta_g, k1, rho1, cp1, runTime.value()).value()
			 << "\t" 
			 << mag(radius.value() - bubbleRadius(beta_g, k1, rho1, cp1, runTime.value()).value())/
					mag(bubbleRadius(beta_g, k1, rho1, cp1, runTime.value()).value() + VSMALL)*100
	  		 << endl;
		
    }

    Info<< "End" << endl;

    return 0;
}

// ************************************************************************* //
