#include "XenPhysicsList.hh"
#include "G4HadronElasticProcess.hh"
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPInelasticData.hh"
#include "G4NeutronIsotopeProduction.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4HadronFissionProcess.hh"
#include "G4NeutronHPFission.hh"
#include "G4NeutronHPFissionData.hh"
#include "G4ProcessManager.hh"

XenPhysicsList::XenPhysicsList() : G4VUserPhysicsList()
{
  defaultCutValue = 1.0*cm;
  SetVerboseLevel(2);
}
XenPhysicsList::~XenPhysicsList()
{

}
void XenPhysicsList::ConstructProcess()
{
  //AddTransportation();
  //ConstructEM();
  //ConstructGeneral();
//   ConstructOp();   // best place to (de)activate optical photon tracking   (part 1 of 2)
  ConstructHad();
//   ConstructLeptHad();
//   ConstructIonFix();
//   ConstructNeutronFix();

}
void XenPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program.

  //ConstructBosons();
  //ConstructLeptons();
  //ConstructMesons();
  //ConstructBaryons();
  //ConstructIons();
}
void XenPhysicsList::ConstructHad()
{

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "neutron") {

//       if (control->UseHPNeutrons()) {

      G4cout << "High precision neutron models chosen" << G4endl;


      // Elastic process
      G4HadronElasticProcess* el_process = new G4HadronElasticProcess();
      pmanager->AddDiscreteProcess(el_process);

      // High precision model and data below 20 MeV
      G4NeutronHPElastic* hpel_model = new G4NeutronHPElastic();
      G4NeutronHPElasticData* el_data = new G4NeutronHPElasticData();
      hpel_model->SetMaxEnergy(19.9*MeV);
      el_process->AddDataSet(el_data);
      el_process->RegisterMe(hpel_model);







      //*************************************************************************
      // Inelastic process ******************************************************

      G4NeutronInelasticProcess* inel_process = new G4NeutronInelasticProcess();

      // High precision model and data below 20 MeV
      G4NeutronHPInelastic* hpinel_model = new G4NeutronHPInelastic();
      G4NeutronHPInelasticData* hpinel_data = new G4NeutronHPInelasticData();

      hpinel_model->SetMaxEnergy(19.9*MeV);

      inel_process->AddDataSet(hpinel_data);
      inel_process->RegisterMe(hpinel_model);

      G4NeutronIsotopeProduction* theNeutronIsotopeProduction = new G4NeutronIsotopeProduction;
      //inel_process->RegisterIsotopeProductionModel(theNeutronIsotopeProduction);//TODO: re add this line
      //G4HadronicProcess::EnableIsotopeProductionGlobally();

      pmanager->AddDiscreteProcess(inel_process);
      //*************************************************************************



      //Capture process
      G4HadronCaptureProcess* cap_process = new G4HadronCaptureProcess();

      // High precision model and data below 20 MeV
      G4NeutronHPCapture* hpcap_model = new G4NeutronHPCapture();
      G4NeutronHPCaptureData* hpcap_data = new G4NeutronHPCaptureData();
      hpcap_model->SetMaxEnergy(19.9*MeV);
      cap_process->AddDataSet(hpcap_data);
      cap_process->RegisterMe(hpcap_model);

//
      pmanager->AddDiscreteProcess(cap_process);


      // Fission process
      G4HadronFissionProcess* fis_process = new G4HadronFissionProcess();

      // High precision model and data below 20 MeV
      G4NeutronHPFission* hpfis_model = new G4NeutronHPFission();
      G4NeutronHPFissionData* hpfis_data = new G4NeutronHPFissionData();
      hpel_model->SetMaxEnergy(19.9*MeV);
      fis_process->AddDataSet(hpfis_data);
      fis_process->RegisterMe(hpfis_model);

      pmanager->AddDiscreteProcess(fis_process);

    }
  }  // while
}
