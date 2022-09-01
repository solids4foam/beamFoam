if (true)
{
    surfaceVectorField DThetaf = fvc::interpolate(DTheta);  

    surfaceScalarField magDThetaf = mag(DThetaf) + SMALL;

    surfaceTensorField DThetaHat = spinTensor(DThetaf);

    DLambda =
        II + (Foam::sin(magDThetaf)/magDThetaf)*DThetaHat
      + ((1.0-Foam::cos(magDThetaf))/sqr(magDThetaf))*(DThetaHat & DThetaHat);
}
