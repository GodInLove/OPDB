package edu.scu.service;

import edu.scu.dao.OperonMapper;
import edu.scu.entity.Operon;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.io.File;
import java.io.IOException;
import java.util.List;

@Service
public class OperonServiceImpl implements OperonService{
    @Autowired
    private OperonMapper operonMapper;

    @Override
    public String isNeed(String srr_num){
        Operon operon = operonMapper.getSrrnum(srr_num);
        if (operon == null){
            return "need";
        }
        return  "no_need";
    }

    @Override
    public String predictOP(String srr_num, String kegg_id , String method)  {
        Integer id_method = 0;
        if (method.equals("CONDOP")){
            id_method = 1;
        }
        String path = "/home/lyd/Documents/OPDB/PyOpdb/";
        String cmd = "python3.5 " + path + "main.py -i " + srr_num +" -o " + path + "test_results -m " + id_method.toString() + " -k " + kegg_id + " > " + path + srr_num + ".log 2>&1";
        File dir = new File(path);
        System.out.println(cmd);
        try {
            Process process;
            process = Runtime.getRuntime().exec(cmd,null,dir);
            int status = process.waitFor();
            if (status != 0){
                System.err.println("Failed to call shell's command");
                return "wrong";
            }
        } catch (IOException e) {
            return e.getMessage() + "\n" + cmd;
        } catch (InterruptedException e) {
            return e.getMessage();
        }
        operonMapper.addOperon(kegg_id,srr_num);
        return "done";
    }
}
