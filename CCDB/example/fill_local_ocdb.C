void fill_local_ocdb()
{
  o2::ccdb::Manager* cdb = o2::ccdb::Manager::Instance();
  cdb->setDefaultStorage("local://O2CDB");
  TH1F* h = 0;
  for (int run = 2000; run < 2100; run++) {
    cdb->setRun(run);
    h = new TH1F(Form("histo for %d run", run), "gauss", 100, -5, 5);
    h->FillRandom("gaus", 1000);
    o2::ccdb::ConditionId* id = new o2::ccdb::ConditionId("DET/Calib/Histo", run, run, 1, 0);
    o2::ccdb::ConditionMetaData* md = new o2::ccdb::ConditionMetaData();
    cdb->putObject(h, *id, md);
    h = 0;
  }
}
