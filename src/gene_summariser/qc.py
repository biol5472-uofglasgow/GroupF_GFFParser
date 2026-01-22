from gene_summariser.models import Transcript


#QC for phase incosistency
def check_cds_phase_errors(transcript:Transcript) -> list[str]:
    flags:list[str]=[]
    cds_list=transcript.cds_features

    if len(cds_list)<2:
        return flags
    
    for i in range (len(cds_list)-1):
        current_cds= cds_list[i]
        next_cds=cds_list[i+1]

        cds_length= current_cds.length

        expected_next_phase=(current_cds.phase+ cds_length)%3
        if next_cds.phase != expected_next_phase:
            flags.append("CDS_Phase_Incosistent")
            break
    
    return flags
