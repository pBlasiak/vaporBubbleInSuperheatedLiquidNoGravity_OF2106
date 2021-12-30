IOdictionary vaporBubbleInSuperheatedLiquidNoGravityDict
(
	IOobject
	(
		"vaporBubbleInSuperheatedLiquidNoGravityDict",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ_IF_MODIFIED,
		IOobject::NO_WRITE
	)
);

// The integrals are solved using rectangular method
// Number of rectangulars
const label N(vaporBubbleInSuperheatedLiquidNoGravityDict.lookupOrDefault("N", 200000));

// Automatic set of starting value for beta_g
dimensionedScalar beta_g("beta_g", dimless, vaporBubbleInSuperheatedLiquidNoGravityDict.lookupOrDefault("beta_g", 0.1));

// Free stream temperature
const dimensionedScalar Tinf("Tinf", dimTemperature, vaporBubbleInSuperheatedLiquidNoGravityDict);

// Bubble radius
const dimensionedScalar R("R", dimLength, vaporBubbleInSuperheatedLiquidNoGravityDict);

// Radius up to which calculations will be performed
const dimensionedScalar Rend("Rend", dimLength, vaporBubbleInSuperheatedLiquidNoGravityDict);

// Radius origin
const dimensionedVector Rorigin("Rorigin", dimLength, vaporBubbleInSuperheatedLiquidNoGravityDict.lookupOrDefault("Rorigin", vector(0,0,0)));

// The number of radius segments from R to Rend
const label Rsegments(vaporBubbleInSuperheatedLiquidNoGravityDict.lookupOrDefault("Rsegments", 100));

Switch calcTiniAnalytical(vaporBubbleInSuperheatedLiquidNoGravityDict.lookupOrDefault("calcTiniAnalytical", false));

// The direction of gradAlpha calculation
const label gradAlphaCalcDir(vaporBubbleInSuperheatedLiquidNoGravityDict.lookupOrDefault("gradAlphaCalcDir", 0));
