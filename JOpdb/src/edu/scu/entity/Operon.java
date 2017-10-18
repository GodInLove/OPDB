package edu.scu.entity;

/**
 * Created by zhaoliang on 2017/10/13.
 */
public class Operon {
    private Integer Id;
    private String kegg_id;
    private String srr_num;
    private String fna_path;
    private String gff_path;
    private String bam_path;
    private String txt_path;
    private String wig_path;

    public Operon(){

    }

    public Operon(String kegg_id, String srr_num, String fna_path, String gff_path, String bam_path, String txt_path, String wig_path) {
        this.kegg_id = kegg_id;
        this.srr_num = srr_num;
        this.fna_path = fna_path;
        this.gff_path = gff_path;
        this.bam_path = bam_path;
        this.txt_path = txt_path;
        this.wig_path = wig_path;
    }

    public Integer getId() {
        return Id;
    }

    public void setId(Integer id) {
        Id = id;
    }

    public String getKegg_id() {
        return kegg_id;
    }

    public void setKegg_id(String kegg_id) {
        this.kegg_id = kegg_id;
    }

    public String getSrr_num() {
        return srr_num;
    }

    public void setSrr_num(String srr_num) {
        this.srr_num = srr_num;
    }

    public String getFna_path() {
        return fna_path;
    }

    public void setFna_path(String fna_path) {
        this.fna_path = fna_path;
    }

    public String getGff_path() {
        return gff_path;
    }

    public void setGff_path(String gff_path) {
        this.gff_path = gff_path;
    }

    public String getBam_path() {
        return bam_path;
    }

    public void setBam_path(String bam_path) {
        this.bam_path = bam_path;
    }

    public String getTxt_path() {
        return txt_path;
    }

    public void setTxt_path(String txt_path) {
        this.txt_path = txt_path;
    }

    public String getWig_path() {
        return wig_path;
    }

    public void setWig_path(String wig_path) {
        this.wig_path = wig_path;
    }

}
