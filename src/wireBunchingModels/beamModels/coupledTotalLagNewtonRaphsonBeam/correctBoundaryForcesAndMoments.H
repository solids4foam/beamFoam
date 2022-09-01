
// Forces
{
    forAll(W_.boundaryField(), patchI)
    {
        if
        (
            isA<forceBeamDisplacementNRFvPatchVectorField>
            (
                W_.boundaryField()[patchI]
            )
        )
        {
            const forceBeamDisplacementNRFvPatchVectorField& pW =
                refCast<const forceBeamDisplacementNRFvPatchVectorField>
                (
                    W_.boundaryField()[patchI]
                );

            explicitQ_.boundaryField()[patchI] = pW.force();
            
            // Info << "force: " << pW.force() << endl;
            // Info << "force: " << explicitQ_.boundaryField()[patchI] << endl;
        }
    }
}

// Moments
{
    forAll(Theta_.boundaryField(), patchI)
    {
        if
        (
            isA<momentBeamRotationNRFvPatchVectorField>
            (
                Theta_.boundaryField()[patchI]
            )
        )
        {
            const momentBeamRotationNRFvPatchVectorField& pTheta =
                refCast<const momentBeamRotationNRFvPatchVectorField>
                (
                    Theta_.boundaryField()[patchI]
                );

            explicitM_.boundaryField()[patchI] = pTheta.moment();

            // Info << "moment: " << pTheta.moment() << endl;
        }
    }
}
