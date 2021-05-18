/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dynamicSmagorinsky.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void dynamicSmagorinsky<BasicTurbulenceModel>::correctNut
(
    const tmp<volTensorField>& gradU
)
{
    this->nut_ = max(Cs_*sqr(this->delta())*mag(dev(symm(gradU))),-1.0*this->nu());
    //this->nut_ = (0.0*Cs_*sqr(this->delta())*mag(dev(symm(gradU))));
    this->nut_.correctBoundaryConditions();
    
    
    fv::options::New(this->mesh_).correct(this->nut_);
    
    BasicTurbulenceModel::correctNut();
}
template<class BasicTurbulenceModel>
void dynamicSmagorinsky<BasicTurbulenceModel>::correctPrtt
(
    const tmp<volTensorField>& gradU
)
{    


  this->Prtt_ =  max(min((Cs_/Cprt_),1.2), 0.6);
 
    this->Prtt_.correctBoundaryConditions();
   //Info<<Prtt_<<endl;
    fv::options::New(this->mesh_).correct(this->Prtt_);
    BasicTurbulenceModel::correctPrtt();
}

template<class BasicTurbulenceModel>
void dynamicSmagorinsky<BasicTurbulenceModel>::correctNut()
{
    correctNut(fvc::grad(this->U_));
}
template<class BasicTurbulenceModel>
void dynamicSmagorinsky<BasicTurbulenceModel>::correctPrtt()
{
    correctPrtt(fvc::grad(this->U_));
}
void correctPrtt(const tmp<volTensorField>& gradU);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
dynamicSmagorinsky<BasicTurbulenceModel>::dynamicSmagorinsky
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const volScalarField& T,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    LESeddyViscosity<BasicTurbulenceModel>
    (
        type,
        alpha,
        rho,
        U,
        T,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", this->U_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    Cs_
    (
        IOobject
        (
            IOobject::groupName("Cs", this->U_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar ("Cs", dimless,SMALL)
    ),
    
   
    
    Cprt_
    (
        IOobject
        (
            IOobject::groupName("Cprt", this->U_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar ("Cprt", dimless,SMALL)
    ),
 
   
    Prtt_
    (
        IOobject
        (
            IOobject::groupName("Prtt", this->U_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
       dimensionedScalar ("Prtt", dimless,SMALL)
    ),
    
    Prtt1_
    (
        IOobject
        (
            IOobject::groupName("Prtt1", this->U_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
       dimensionedScalar ("Prtt1", dimless,SMALL)
    ),
    
    
   

   
    
    simpleFilter_(this->mesh_),
    filterPtr_(LESfilter::New(this->mesh_, this->coeffDict())),
    filter_(filterPtr_())
    
{
   // bound(k_, this->kMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool dynamicSmagorinsky<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        filter_.read(this->coeffDict());        

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void dynamicSmagorinsky<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const surfaceScalarField& phi = this->phi_;
    const volVectorField& U = this->U_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    LESeddyViscosity<BasicTurbulenceModel>::correct();

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    volSymmTensorField S(dev(symm(gradU)));
    volScalarField magS(mag(S));

    volVectorField Uf(filter_(U));

    volSymmTensorField Sf(filter_(S));  

    
    volScalarField magSf(mag(Sf));
          
    volSymmTensorField LL
    (
        simpleFilter_(dev(filter_(sqr(U)) - (sqr(filter_(U)))))
    );
   

    volSymmTensorField MM
    (
        simpleFilter_(sqr(this->delta())*(filter_(magS*S) - 4.0*magSf*Sf))
    );
    
    

    Cs_= (
        simpleFilter_((LL && MM))
       /(
            simpleFilter_(magSqr(MM))
          + dimensionedScalar("small", sqr(MM.dimensions()), VSMALL)
        )
    );

  
    volScalarField KK =
    (filter_(magSqr(U)) - magSqr(filter_(U)));
    
    volScalarField mm
    (
        sqr(this->delta())*(sqr(mag(Sf)) - filter_(sqr(magS)))
       
    );

    volScalarField mmmm = fvc::average(magSqr(mm));
    mmmm.max(VSMALL);

    k_ = fvc::average(KK*mm)/mmmm * sqr(this->delta())*magSqr(S);

    correctNut(gradU);
    

   const volScalarField& T = this->T_;
    tmp<volVectorField> tgradT(fvc::grad(T));
    
    const volVectorField& gradT = tgradT();
  
   
     volVectorField PP 
     
     (
    simpleFilter_(filter_((T)*(U)) - (filter_(U)*filter_(T))));
    
 
    volVectorField RR =
    (
       simpleFilter_(sqr(this->delta())*( (4.0*magSf*filter_(gradT))- filter_((magS)*(gradT
    )))));
    
  
    Cprt_=
    (
        simpleFilter_(0.5*(PP && RR))
       /(
            simpleFilter_(magSqr(RR)) + dimensionedScalar("small", sqr(RR.dimensions()), SMALL)
          )
        
    );
    Prtt1_=(Cs_/Cprt_);
    //Info<<Cprt_<<endl;
    
    correctPrtt(gradU);
    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
