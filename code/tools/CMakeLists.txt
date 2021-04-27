set (TOOLS
  bm_bnlb_tool
  bm_extract_tool
  bm_stitch_tool
)

foreach (aux ${TOOLS})
 add_executable (${aux} ${aux}.cpp)
 target_link_libraries(${aux} binmat ${EXTLIB}) 
endforeach (aux)
