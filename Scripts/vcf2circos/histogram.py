from pprint import pprint
from typing import Generator
from vcf2circos.plotcategories.plotconfig import Plotconfig
from vcf2circos.utils import timeit, generate_hovertext_var, chr_valid

from os.path import join as osj
import os
import pandas as pd
from itertools import chain, repeat
from collections import OrderedDict, Counter


class Histogram_(Plotconfig):
    """
    It need to create one histogram for each SV event FOR EACH SV height (from 0 copy number to 5), which will create the grey band between color dor
    """

    def __init__(self, plotconfig):
        self.plotconfig = plotconfig
        self.variants_position = self.config_ring["position"]
        self.variants_ring_space = self.config_ring["space"]
        self.variants_ring_height = self.config_ring["height"]
        # corresponding to SNV InDel height 7th ring (after 0 to 5 copy number height)
        self.radius = {
            "R0": self.variants_position
            + (max(self.rangescale) * self.variants_ring_space)
            + ((max(self.rangescale) + 1) * self.variants_ring_height),
            "R1": self.variants_position
            + (max(self.rangescale) * self.variants_ring_space)
            + ((max(self.rangescale) + 2) * self.variants_ring_height),
        }
        # print("#Range", self.rangescale)
        self.hovertextformat = (
            ' "<b>{}:{}-{}</b><br>{}<br><br>{}".format(a[i,0], a[i,1], a[i,2], a[i,6], a[i,8])'
        )

        # self.hovertextformat = ""
        self.trace = {
            "hoverinfo": "text",
            "mode": "markers",
            "marker": {
                "size": 5,
                "symbol": 0,
                "color": self.colors["INTERMEDIATE"],
                "opacity": 0.1,
            },
        }
        self.layout = {
            "type": "path",
            "opacity": 1,
            "fillcolor": self.colors["INTERMEDIATE"],
            "line": {"color": self.colors["INTERMEDIATE"], "width": 5},
        }
        # TODO tile same as cytobandinfo in vfreader
        self.cytoband_conf = pd.read_csv(
            osj(
                self.options["Static"],
                "Assembly",
                self.options["Assembly"],
                "cytoband_" + self.options["Assembly"] + "_chr_infos.txt.gz",
            ),
            sep="\t",
            header=0,
            compression="infer",
        )
        self.cytoband_data = {
            "show": "True",
            "file": {
                "path": "",
                "header": "infer",
                "sep": "\t",
                "dataframe": {"orient": "columns", "data": None},
            },
            "colorcolumn": 4,
            "radius": {"R0": 1, "R1": 1.1},
            "hovertextformat": ' "<b>{}:{}-{}<br>{}</b>".format(a[i,0], a[i,1], a[i,2], a[i,5])',
            "trace": {
                "uid": "cytoband_tile",
                "hoverinfo": "text",
                "mode": "markers",
                "marker": {
                    "size": 0,
                    "symbol": 0,
                    "color": None,
                    "opacity": 0,
                },  # 8
            },
            "layout": {
                "type": "path",
                "layer": "above",
                "opacity": 0,
                "line": {"color": None, "width": 0},
            },
        }
        self.only_overlapping_ = self.options["Genes"]["only_snv_in_sv_genes"]

    def __getattr__(self, item):
        if hasattr(self.plotconfig, item):
            return getattr(self.plotconfig, item)

    def dict_to_str(self, info_field: list) -> Generator:
        for info_dict in info_field:
            yield ";".join([str(key) + "=" + str(value) for key, value in info_dict.items()])

    def only_snv_indels_in_sv(self, data):
        df_data = pd.DataFrame.from_dict(data).astype(
            {
                "chr_name": str,
                "start": int,
                "end": int,
                "val": int,
                "ref": str,
                "alt": str,
                "type": str,
                "color": str,
                "hovertext": str,
                "symbol": int,
                "genes": str,
                "exons": str,
            }
        )
        snv_indels_df = df_data[df_data["type"].isin(["SNV", "INDEL", "OTHER"])]
        sv_df = df_data[~df_data["type"].isin(["SNV", "INDEL", "OTHER"])]
        # pd.set_option("display.max_columns", None)
        # pd.set_option("display.width", None)
        # pd.set_option("display.max_colwidth", -1)
        # SNV / INDEL
        snv_indel_not_overlapp = []
        snv_indel_overlapp = []
        for i, si in snv_indels_df.iterrows():
            # WHOLE
            for j, v in sv_df.iterrows():
                if (
                    si["chr_name"] == v["chr_name"]
                    and si["start"] > v["start"]
                    and si["start"] < v["end"]
                ):
                    snv_indel_overlapp.append([si["chr_name"], si["start"], si["ref"], si["alt"]])
                # else:
                #    if (
                #        not [si["chr_name"], si["start"], si["end"]]
                #        in snv_indel_not_overlapp
                #    ):
                #        snv_indel_not_overlapp.append(
                #            [si["chr_name"], si["start"], si["end"]]
                #        )
        # print(snv_indel_not_overlapp)
        return snv_indel_overlapp

    def adapt_data(self, cn: int) -> dict:
        d_file = {
            "path": "",
            "header": "infer",
            "sep": "\t",
            "dataframe": {"orient": "columns", "data": None},
        }
        data = {
            "chr_name": [],
            "start": [],
            "end": [],
            "val": [],
            "ref": [],
            "alt": [],
            "type": [],
            "color": [],
            "hovertext": [],
            "symbol": [],
            "genes": [],
            "exons": [],
        }

        df_data = self.df_data.loc[self.df_data["CopyNumber"] == cn]
        start = []
        stop = []
        ref = []
        alt = []
        for items in list(
            self.extract_start_stop_ref_alt(
                df_data["Record"].to_list(),
                df_data["Variants"].to_list(),
                df_data["Variants_type"].to_list(),
            )
        ):
            try:
                # DEBUGG
                # print(*items)
                # for val in items:
                #    if isinstance(val, list):
                #        print(type(val[0]))
                #    else:
                #        print(type(val))
                # exit()
                start.append(items[0])
                stop.append(items[1])
                ref.append(items[2])
                alt.append(str(items[3][0]))
            except IndexError:
                print("ERROR ", items)
                exit()

        data["chr_name"].extend(df_data["Chromosomes"].to_list())
        data["start"].extend(start)
        data["end"].extend(stop)
        data["val"].extend(list(repeat(2, len(df_data.index))))
        data["ref"].extend(ref)
        data["alt"].extend(alt)
        data["type"].extend(df_data["Variants_type"].to_list())
        data["color"].extend(list(repeat(self.colors["INTERMEDIATE"], len(df_data.index))))
        # data["hovertext"].extend(list(itertools.repeat("", len(df_data.index))))
        data["hovertext"].extend(
            list(
                generate_hovertext_var(
                    df_data["Variants"],
                    full_annot=20,
                    true_annot=self.options["Variants"]["annotations"]["fields"],
                )
            )
        )
        # print(data["hovertext"])
        # data["hovertext"].extend(
        #    [
        #        "Genes ("
        #        + str(len(record.split(",")))
        #        + "): "
        #        + ",".join(record.split(",")[:5])
        #        for record in df_data["Genes"].to_list()
        #    ]
        # )
        data["symbol"].extend(list(repeat(0, len(df_data.index))))
        data["genes"].extend(df_data["Genes"].to_list())
        data["exons"].extend(list(repeat("", len(df_data.index))))
        # data["info"].extend(list(self.dict_to_str(df_data["Variants"].to_list())))

        d_file["dataframe"]["data"] = data
        return d_file

    def histo_cnv_level(self, cn: int) -> dict:
        d = {}
        d_file = self.adapt_data(cn)
        d["show"] = "True"
        d["customfillcolor"] = "False"
        d["file"] = d_file
        d["sortbycolor"] = "False"
        d["colorcolumn"] = 7
        try:
            radius = (
                self.rangescale[cn]
                + self.rangescale[cn]
                + self.options["Variants"]["rings"]["height"]
            ) / 2
        except TypeError:
            print("ERROR radius cnv level")
            exit()
        d["radius"] = {
            "R0": radius,
            "R1": radius,
        }
        d["hovertextformat"] = self.hovertextformat
        d["trace"] = {
            "hoverinfo": "text",
            "mode": "markers",
            "marker": {
                "size": 5,
                "symbol": d_file["dataframe"]["data"]["symbol"],
                "color": self.colors["INTERMEDIATE"],
                "opacity": 0.3,
            },
            "uid": "cnv_scatter_level_" + str(cn),
        }
        d["layout"] = {
            "type": "path",
            "layer": "above",
            "opacity": 0.3,
            "fillcolor": "red",
            "line": {"color": self.colors["INTERMEDIATE"], "width": 5},
        }
        return d

    def cytoband_tile(self, cytoband_data):
        """
        Needed to have annotations in cytoband ring, invisible histogram in background of each band
        """
        dico_cyto = self.cytoband_data.copy()
        cyto = {}
        cyto["chr_name"] = cytoband_data["chr_name"]
        cyto["start"] = cytoband_data["start"]
        cyto["end"] = cytoband_data["end"]
        # Remember to have val column in data otherwise it leads to crash]
        cyto["val"] = list(repeat(1, len(cytoband_data["chr_name"])))
        cyto["band_color"] = list(repeat(self.colors["CYTOBAND"], len(cytoband_data["chr_name"])))
        cyto["band"] = cytoband_data["band"]
        # Cytoband tiles 3  need fill data
        dico_cyto["file"]["dataframe"]["data"] = cyto

        dico_cyto["layout"]["line"]["color"] = cyto["band_color"]
        dico_cyto["trace"]["marker"]["color"] = cyto["band_color"]
        # print(dico_cyto)
        # exit()
        return dico_cyto

    def merge_options(self, cytoband_data: dict) -> list:
        """
        func handle math and geometry need to take data in specific order
        chr start stop val OTHERWIS TROUBLE
        """
        # exit()

        whole_cn = []
        # Histo_cnv_level
        for cn in list(set(self.data["CopyNumber"])):
            res = self.histo_cnv_level(cn)
            whole_cn.append(res)

        if self.only_overlapping_:
            print("#[INFO] SNV / indels Overlapping SV only")
            whole_var = {
                "chr_name": [],
                "start": [],
                "end": [],
                "val": [],
                "ref": [],
                "alt": [],
                "type": [],
                "color": [],
                "hovertext": [],
                "symbol": [],
                "genes": [],
                "exons": [],
            }
            for dic in whole_cn:
                for key, val in dic["file"]["dataframe"]["data"].items():
                    whole_var[key].extend(val)
            snv_indel_overlapp = self.only_snv_indels_in_sv(whole_var)
            # print(snv_indel_overlapp)
            for dico in whole_cn:
                # print(dico["trace"]["uid"])
                if dico["trace"]["uid"] == "cnv_scatter_level_6":
                    df_ = pd.DataFrame.from_dict(dico["file"]["dataframe"]["data"])
                    # if at least one snv indel overlap a sv
                    if snv_indel_overlapp:
                        for wr in snv_indel_overlapp:
                            df_ = df_.loc[
                                (df_["chr_name"] == wr[0])
                                & (df_["start"] == wr[1])
                                & (df_["ref"] == wr[2])
                                & (df_["alt"] == wr[3])
                            ]
                            dico["file"]["dataframe"]["data"] = df_.to_dict("list")
                    else:
                        dico["file"]["dataframe"]["data"] = {
                            "chr_name": [],
                            "start": [],
                            "end": [],
                            "val": [],
                            "ref": [],
                            "alt": [],
                            "type": [],
                            "color": [],
                            "hovertext": [],
                            "symbol": [],
                            "genes": [],
                            "exons": [],
                        }
                # remaining_var.append(self.remove_snv_(snv_indel_overlapp))
        # Extra
        whole_cn.extend(self.generate_extra_plots_from_df())

        # Genes plots
        whole_cn.append(self.histo_genes())

        # cytoband tiles
        whole_cn.append(self.cytoband_tile(cytoband_data))
        return whole_cn

    def remove_snv_(self, list_to_remove):
        index_keep = []
        index_to_rm = []
        for vars in list_to_remove:
            # print(vars)
            for i, record in enumerate(self.data["Record"]):
                if (
                    record.CHROM == vars[0]
                    and record.POS == vars[1]
                    and record.REF == vars[2]
                    and str(record.ALT[0]) == vars[3]
                ):
                    index_keep.append(i)
        for j, rd in enumerate(self.data["Record"]):
            # print(rd.var_type)
            # print(rd)
            # print("\n")
            if rd.var_type == "snp" or rd.var_type == "indel":
                # print(rd)
                # print(j)
                if j not in index_keep:
                    # print("index to remove: " + str(j))
                    # print(record)
                    index_to_rm.append(j)

        self.df_data.drop(index=index_to_rm, inplace=True)
        self.data = self.df_data.to_dict("list")
        return self.data

    def process_gene_list(self, genes_list: list) -> Generator:
        for record in genes_list:
            if record:
                yield record.split(",")

    def histo_genes(self) -> dict:
        data = {}
        dico = {}
        # remove empty gene, df_data attribute of class basic data from plot config Parents class
        # gene_list = list(filter(lambda x: x != "", self.df_data["Genes"]))
        gene_list = list(
            set(
                list(
                    map(
                        str,
                        chain.from_iterable(list(self.process_gene_list(self.df_data["Genes"]))),
                    )
                )
            )
        )
        # print(*self.df_genes.columns)
        # print(self.df_genes.head())
        ## select genes in or batch of variations (from refeseq assembly)
        df_filter = self.df_genes.loc[
            (self.df_genes["gene"].isin(gene_list)) & (self.df_genes["chr_name"]).isin(chr_valid())
        ]

        # print(*gene_list)
        # print(self.df_data["Genes"].head())

        # Set color
        for fields in df_filter.columns:
            if fields != "transcript":
                if fields == "color":
                    # data[fields] = list(self.morbid_genes(df_filter["gene"]))
                    data[fields] = list(repeat(self.colors["GENES"], len(df_filter.index)))
                else:
                    data[fields] = df_filter[fields].to_list()
        # pprint(data, sort_dicts=False)
        dico["file"] = {
            "path": "",
            "header": "infer",
            "sep": "\t",
            "dataframe": {"orient": "columns", "data": data},
        }
        dico["show"] = self.show
        dico["colorcolumn"] = 4
        dico["radius"] = {"R0": 0.96, "R1": 0.96}
        dico[
            "hovertextformat"
        ] = ' "<b>{}:{}-{}<br>Gene: {}</b><br>".format(a[i,0], a[i,1], a[i,2], a[i,5])'
        dico["trace"] = {
            "uid": "genes",
            "hoverinfo": "text",
            "mode": "markers",
            "marker": {
                "size": 3,
                "symbol": 0,
                "color": data["color"],
                "opacity": 1,
            },
        }
        dico["layout"] = {
            "type": "path",
            "layer": "above",
            "opacity": 0.2,
            "line": {"color": data["color"], "width": 3},
        }
        return dico

    def genes_omim_morbid(self):
        """ "
        If it's a morbid gene it will be colored in red in circos gene level
        done in static file in genes.<assembly>
        """
        pass

    def extract_start_stop_ref_alt(
        self, record: list, info_field: list, variant_type: list
    ) -> Generator:
        """
        修复版本 - 处理特殊的VCF INFO格式 (INFO=END=2652177|SVLEN=-2600|SVTYPE=CNV)
        """
        for i, info_dict in enumerate(info_field):
            if variant_type[i] not in ["OTHER", "SNV", "INDEL"]:
                
                # 步骤1: 预处理INFO字段，处理特殊格式
                print(f"处理记录 {i+1}: {record[i].CHROM}:{record[i].POS}")
                print(f"变异类型: {variant_type[i]}")
                print(f"原始INFO: {dict(record[i].INFO)}")
                
                # 检查每个INFO字段的值，寻找特殊格式
                processed_info = {}
                
                for key, value in record[i].INFO.items():
                    if isinstance(value, list):
                        # 处理列表格式的INFO字段
                        print(f"检测到列表格式的INFO字段 {key}: {value}")
                        
                        # 如果key是'INFO'且value是包含'='的字符串列表
                        if key == 'INFO' and all('=' in str(item) for item in value):
                            print("检测到特殊的INFO列表格式，开始解析...")
                            
                            for item in value:
                                item_str = str(item)
                                if '=' in item_str:
                                    try:
                                        field_name, field_value = item_str.split('=', 1)
                                        
                                        # 转换数值
                                        if field_value.lstrip('-+').replace('.', '').isdigit():
                                            processed_info[field_name] = int(float(field_value))
                                        else:
                                            processed_info[field_name] = field_value
                                            
                                        print(f"  提取字段: {field_name} = {processed_info[field_name]}")
                                    except ValueError as e:
                                        print(f"  解析错误: {item_str} -> {e}")
                        else:
                            # 处理其他列表格式
                            str_val = str(value[0]) if len(value) > 0 else ""
                            
                            # 检测包含管道符和等号的格式
                            if '|' in str_val and '=' in str_val:
                                print(f"检测到特殊格式在字段 {key}: {str_val}")
                                
                                # 解析特殊格式，分割并重新解析各个字段
                                parts = str_val.split('|')
                                for part in parts:
                                    if '=' in part:
                                        try:
                                            field_name, field_value = part.split('=', 1)
                                            
                                            # 转换数值
                                            if field_value.lstrip('-+').replace('.', '').isdigit():
                                                processed_info[field_name] = int(float(field_value))
                                            else:
                                                processed_info[field_name] = field_value
                                                
                                            print(f"  提取字段: {field_name} = {processed_info[field_name]}")
                                        except ValueError as e:
                                            print(f"  解析错误: {part} -> {e}")
                            
                            # 处理包含管道符的单一值（如 "100|200"）
                            elif '|' in str_val and '=' not in str_val:
                                print(f"检测到管道符分隔的值: {key} = {str_val}")
                                try:
                                    # 取第一个值
                                    first_value = str_val.split('|')[0]
                                    if first_value.lstrip('-+').replace('.', '').isdigit():
                                        processed_info[key] = int(float(first_value))
                                    else:
                                        processed_info[key] = first_value
                                    print(f"  使用第一个值: {key} = {processed_info[key]}")
                                except (ValueError, IndexError):
                                    print(f"  无法解析管道符分隔的值: {str_val}")
                            else:
                                # 标准列表处理
                                processed_info[key] = value
                    
                    elif isinstance(value, str):
                        # 处理字符串格式
                        if '|' in value and '=' in value:
                            print(f"检测到特殊格式在字段 {key}: {value}")
                            
                            # 解析特殊格式，分割并重新解析各个字段
                            parts = value.split('|')
                            for part in parts:
                                if '=' in part:
                                    try:
                                        field_name, field_value = part.split('=', 1)
                                        
                                        # 转换数值
                                        if field_value.lstrip('-+').replace('.', '').isdigit():
                                            processed_info[field_name] = int(float(field_value))
                                        else:
                                            processed_info[field_name] = field_value
                                            
                                        print(f"  提取字段: {field_name} = {processed_info[field_name]}")
                                    except ValueError as e:
                                        print(f"  解析错误: {part} -> {e}")
                        
                        # 处理包含管道符的单一值（如 "100|200"）
                        elif '|' in value and '=' not in value:
                            print(f"检测到管道符分隔的值: {key} = {value}")
                            try:
                                # 取第一个值
                                first_value = value.split('|')[0]
                                if first_value.lstrip('-+').replace('.', '').isdigit():
                                    processed_info[key] = int(float(first_value))
                                else:
                                    processed_info[key] = first_value
                                print(f"  使用第一个值: {key} = {processed_info[key]}")
                            except (ValueError, IndexError):
                                print(f"  无法解析管道符分隔的值: {value}")
                        else:
                            # 标准字符串
                            processed_info[key] = value
                    else:
                        # 其他类型（int等）
                        processed_info[key] = value
                
                # 步骤2: 更新record的INFO字段
                record[i].INFO.update(processed_info)
                print(f"处理后的INFO: {dict(record[i].INFO)}")
                
                # 步骤3: 标准化剩余的INFO字段中的数值
                for values in ["SV_start", "SV_end", "SVLEN", "END"]:
                    if values in record[i].INFO:
                        current_value = record[i].INFO[values]
                        
                        if isinstance(current_value, list):
                            if len(current_value) > 0:
                                str_val = str(current_value[0])
                                if "|" in str_val:
                                    record[i].INFO[values] = int(str_val.split("|")[0])
                                else:
                                    record[i].INFO[values] = int(str_val)
                        elif isinstance(current_value, str):
                            if "|" in current_value:
                                record[i].INFO[values] = int(current_value.split("|")[0])
                            else:
                                record[i].INFO[values] = int(current_value)
                        # 如果已经是int，保持不变
                        
                        print(f"标准化 {values}: {record[i].INFO[values]}")
                
                # 步骤4: 确定起始和结束位置
                start_pos = None
                end_pos = None
                
                if "SV_start" in record[i].INFO and "SV_end" in record[i].INFO:
                    start_pos = record[i].INFO["SV_start"]
                    end_pos = record[i].INFO["SV_end"]
                    print(f"使用 SV_start/SV_end: {start_pos} -> {end_pos}")
                    
                elif "END" in record[i].INFO:
                    start_pos = int(record[i].POS)
                    end_pos = record[i].INFO["END"]
                    print(f"使用 POS/END: {start_pos} -> {end_pos}")
                    
                elif "SVLEN" in record[i].INFO:
                    start_pos = int(record[i].POS)
                    svlen = record[i].INFO["SVLEN"]
                    end_pos = start_pos + abs(svlen)
                    print(f"使用 SVLEN: {start_pos} + {abs(svlen)} = {end_pos}")
                    
                else:
                    print("警告: 缺少长度信息，尝试推断...")
                    start_pos = int(record[i].POS)
                    
                    # 从REF/ALT推断
                    try:
                        ref_len = len(record[i].REF) if record[i].REF != '.' else 1
                        alt_len = max([len(str(alt)) for alt in record[i].ALT]) if record[i].ALT else 1
                        estimated_len = max(abs(ref_len - alt_len), 100)
                        end_pos = start_pos + estimated_len
                        
                        print(f"  REF长度: {ref_len}")
                        print(f"  ALT长度: {alt_len}")
                        print(f"  推断长度: {estimated_len}")
                        print(f"  推断位置: {start_pos} -> {end_pos}")
                    except Exception as e:
                        print(f"  推断失败: {e}")
                        # 使用默认值
                        end_pos = start_pos + 1000
                        print(f"  使用默认长度: {start_pos} -> {end_pos}")
                
                if start_pos is not None and end_pos is not None:
                    yield (start_pos, end_pos, record[i].REF, record[i].ALT)
                else:
                    print("错误: 无法确定位置信息")
                    yield (int(record[i].POS), int(record[i].POS) + 1000, record[i].REF, record[i].ALT)
            
            # SNV/INDEL处理
            else:
                try:
                    alternate = int(str(max([len(alt) for alt in record[i].ALT])))
                    start_pos = int(str(record[i].POS))
                    end_pos = start_pos + alternate
                    print(f"SNV/INDEL: {start_pos} -> {end_pos}")
                    yield (start_pos, end_pos, record[i].REF, record[i].ALT)
                except Exception as e:
                    print(f"SNV/INDEL处理错误: {e}")
                    yield (int(record[i].POS), int(record[i].POS) + 1, record[i].REF, record[i].ALT)

    def generate_extra_plots_from_df(self):
        extras = []
        if "gc" in self.options["Extra"]:
            # self.gcplus = pd.DataFrame(osj(self.options["static"], #"histogram_pos_chr"))
            gc_pos = osj(
                self.options["Static"],
                "Assembly",
                self.options["Assembly"],
                self.options["Assembly"] + ".histogram_pos_chr.txt",
            )
            gc_neg = osj(
                self.options["Static"],
                "Assembly",
                self.options["Assembly"],
                self.options["Assembly"] + ".histogram_neg_chr.txt",
            )
            # for f in [gc_pos, gc_neg]:
            #    assert os.path.exists(f), (
            #        f + " file does not exists, add in Static folder"
            #    )
            gc_mod = osj(
                self.options["Static"],
                "Assembly",
                self.options["Assembly"],
                self.options["Assembly"] + ".gc5Base.5Mb.notnorm.txt",
            )
            for gc_ in [gc_mod]:
                if os.path.exists(gc_):
                    gc_dict = {
                        "show": "True",
                        "customfillcolor": "False",
                        "file": {
                            "path": gc_,
                            "header": "infer",
                            "sep": "\t",
                        },
                        "sortbycolor": "False",
                        "colorcolumn": "None",
                        "radius": {"R0": 0.88, "R1": 0.92},
                        "hovertextformat": ' "Chromosome: {}<br>Start: {}<br>End: {}    <br>GC mean:{}".format(a[i,0], a[i,1], a[i,2], float(a[i,3])) ',
                        "trace": {
                            "hoverinfo": "text",
                            "mode": "markers",
                            "marker": {"size": 0, "opacity": 0},
                            "uid": "extra_gc",
                        },
                        "layout": {
                            "type": "path",
                            "opacity": 1,
                            "fillcolor": "blue",
                            "line": {"color": "blue", "width": 0},
                        },
                    }
                    extras.append(gc_dict)
                else:
                    print("[WARN GC file don't exist for assembly " + self.options["Assembly"])

        if "mappability_old" in self.options["Extra"]:
            mapp_file = osj(
                self.options["Static"],
                "Assembly",
                self.options["Assembly"],
                self.options["Assembly"] + ".dukeExcludeRegions.csv",
            )
            if not os.path.exists(mapp_file):
                print(
                    "[WARN] Mappability Exclude Region file not in Static folder for assembly "
                    + self.options["Assembly"]
                )
            else:
                data = pd.read_csv(
                    mapp_file,
                    header=0,
                    sep="\t",
                )
                data = data.loc[data["chr_name"].isin(chr_valid())]
                data["val"] = 2
                data["color"] = "red"
                data["ref"] = ""
                data["alt"] = ""
                # data["infos"] = ""
                # data["hovertext"] = ""
                # data["infos_dict"] = ""
                mappa_dict = {
                    "show": "True",
                    "customfillcolor": "False",
                    "file": {
                        "path": "",
                        "header": "infer",
                        "sep": "\t",
                        "dataframe": {
                            "orient": "columns",
                            "data": data.to_dict("list"),
                        },
                    },
                    "sortbycolor": "False",
                    "colorcolumn": 7,
                    "radius": {"R0": 0.80, "R1": 0.84},
                    "hovertextformat": ' "Chromosome: {}<br>Start: {}<br>End: {}<br>Type:   {}".format(a[i,0], a[i,1], a[i,2], a[i,6]) ',
                    "trace": {
                        "hoverinfo": "text",
                        "mode": "markers",
                        "marker": {"size": 0, "opacity": 0},
                        "uid": "extra_mappability",
                    },
                    "layout": {
                        "type": "path",
                        "opacity": 1,
                        "fillcolor": "black",
                        "line": {"color": "black", "width": 0},
                    },
                }
                # print(data)
                extras.append(mappa_dict)
            # print(mappa_dict)
            # exit()
        if "mappability" in self.options["Extra"]:
            filename = [
                bfiles
                for bfiles in os.listdir(
                    osj(self.options["Static"], "Assembly", self.options["Assembly"])
                )
                if bfiles.startswith(self.options["Assembly"] + ".blacklist")
            ][0]
            mapp_file = osj(self.options["Static"], "Assembly", self.options["Assembly"], filename)
            if not os.path.exists(mapp_file):
                print(
                    "[WARN] Blacklist Region file not in Static folder for assembly "
                    + self.options["Assembly"]
                )
            else:
                data = pd.read_csv(mapp_file, header=None, sep="\t", compression="infer")
                data.columns = ["chr_name", "start", "end", "type"]
                data = data.loc[data["chr_name"].isin(chr_valid())]
                data["val"] = 2
                data["color"] = "red"
                data["ref"] = ""
                data["alt"] = ""
                mappa_dict = {
                    "show": "True",
                    "customfillcolor": "False",
                    "file": {
                        "path": "",
                        "header": "infer",
                        "sep": "\t",
                        "dataframe": {
                            "orient": "columns",
                            "data": data.to_dict("list"),
                        },
                    },
                    "sortbycolor": "False",
                    "colorcolumn": 7,
                    "radius": {"R0": 0.80, "R1": 0.84},
                    "hovertextformat": ' "Chromosome: {}<br>Start: {}<br>End: {}<br>Type:   {}".format(a[i,0], a[i,1], a[i,2], a[i,6]) ',
                    "trace": {
                        "hoverinfo": "text",
                        "mode": "markers",
                        "marker": {"size": 0, "opacity": 0},
                        "uid": "extra_mappability",
                    },
                    "layout": {
                        "type": "path",
                        "opacity": 1,
                        "fillcolor": "black",
                        "line": {"color": "black", "width": 0},
                    },
                }
                data = data.loc[
                    :,
                    ["chr_name", "start", "end", "val", "ref", "alt", "type", "color"],
                ]
                # exit()
                # extras.append(mappa_dict)

        if "repeatmasker_old" in self.options["Extra"]:
            assert os.path.exists(
                osj(
                    self.options["Static"],
                    "Assembly",
                    self.options["Assembly"],
                    self.options["Assembly"] + ".repeatmasker.tsv",
                )
            ), (
                "Repeat Masker file "
                + self.options["Assembly"]
                + ".repeatmasker.tsv"
                + " not in Static folder"
            )
            dat = pd.read_csv(
                osj(
                    self.options["Static"],
                    "Assembly",
                    self.options["Assembly"],
                    self.options["Assembly"] + ".repeatmasker.tsv",
                ),
                header=0,
                sep="\t",
                compression="infer",
            )
            data = dat.loc[dat["chr_name"].isin(chr_valid())]
            data["val"] = 2
            data["color"] = self.options["Color"]["MAPPABILITY"]
            data["ref"] = ""
            data["alt"] = ""
            # data["infos"] = ""
            # data["hovertext"] = ""
            # data["infos_dict"] = ""
            data_filter = data.loc[
                :, ["chr_name", "start", "stop", "val", "color", "ref", "alt", "type"]
            ]
            repeat_dict = {
                "show": "True",
                "customfillcolor": "False",
                "file": {
                    "path": "",
                    "header": "infer",
                    "sep": "\t",
                    "dataframe": {
                        "orient": "columns",
                        "data": data_filter.to_dict("list"),
                    },
                },
                "sortbycolor": "False",
                "colorcolumn": 7,
                "radius": {"R0": 0.75, "R1": 0.79},
                "hovertextformat": ' "Chromosome: {}<br>Start: {}<br>End: {}<br>Type:{}".format(a[i,0], a[i,1], a[i,2], a[i,6]) ',
                "trace": {
                    "hoverinfo": "text",
                    "mode": "markers",
                    "marker": {"size": 0, "opacity": 0},
                    "uid": "extra_repeatmasker",
                },
                "layout": {
                    "type": "path",
                    "opacity": 1,
                    "fillcolor": "black",
                    "line": {"color": "black", "width": 0},
                },
            }
            extras.append(repeat_dict)
        if "repeatmasker" in self.options["Extra"]:
            assert os.path.exists(
                osj(
                    self.options["Static"],
                    "Assembly",
                    self.options["Assembly"],
                    self.options["Assembly"] + ".repeatmasker.histo.tsv",
                )
            ), (
                "Repeat Masker file "
                + self.options["Assembly"]
                + ".repeatmasker.histo.tsv"
                + " not in Static folder"
            )
            repeat_dict = {
                "show": "True",
                "customfillcolor": "False",
                "file": {
                    "path": osj(
                        self.options["Static"],
                        "Assembly",
                        self.options["Assembly"],
                        self.options["Assembly"] + ".repeatmasker.histo.tsv",
                    ),
                    "header": "infer",
                    "sep": "\t",
                },
                "sortbycolor": "False",
                "colorcolumn": "",
                "radius": {"R0": 0.75, "R1": 0.79},
                "hovertextformat": ' "Chromosome: {}<br>Start: {}<br>End: {}<br>Type:{}".format(a[i,0], a[i,1], a[i,2])',
                "trace": {
                    "hoverinfo": "text",
                    "mode": "markers",
                    "marker": {"size": 0, "opacity": 0},
                    "uid": "extra_repeatmasker",
                },
                "layout": {
                    "type": "path",
                    "opacity": 1,
                    "fillcolor": "black",
                    "line": {"color": "black", "width": 0},
                },
            }
            extras.append(repeat_dict)
        return extras
