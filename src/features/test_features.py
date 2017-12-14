
import os
import sys

# add the 'src/' directory as one where we can import modules
src_dir = os.path.join('src')
sys.path.append(src_dir)

def test_make_roms_ds():
    
    from features.roms_ds import make_roms_ds
    import os
    
    file_path = os.path.join(os.pardir,'data','raw','waom10_full_forcing','ocean_avg_000[4,5].nc')
    ds = make_roms_ds(file_path)
    
    assert ds.depth.isnull().any() == True, 'depth has no NaN'