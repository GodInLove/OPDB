package edu.scu.dao;

import edu.scu.entity.Operon;
import org.apache.ibatis.annotations.Param;
import org.springframework.stereotype.Repository;

import java.sql.Date;
import java.util.List;

@Repository
public interface OperonMapper {
	public Operon getSrrnum(@Param("srr_num")String srr_num);
	public void addOperon(@Param("kegg_id")String kegg_id, @Param("srr_num")String srr_num);

}
