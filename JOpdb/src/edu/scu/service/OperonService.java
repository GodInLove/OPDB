package edu.scu.service;

import edu.scu.entity.Operon;

import java.io.IOException;
import java.util.List;

public interface OperonService {
    public String isNeed(String srr_num);
    public String predictOP(String srr_num, String kegg_id , String method) throws InterruptedException, IOException;
}
