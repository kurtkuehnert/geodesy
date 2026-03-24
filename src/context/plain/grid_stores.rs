use super::{GLOBALLY_ALLOCATED_GRIDS, GridCollection, Plain, init_grids};
use crate::authoring::*;
use std::{path::PathBuf, sync::Arc};

pub(super) trait GridStore {
    fn get_grid(&self, name: &str, ctx: &Plain) -> Result<Option<Arc<BaseGrid>>, Error>;
}

pub(super) struct MemoryGridStore;
impl GridStore for MemoryGridStore {
    fn get_grid(&self, name: &str, ctx: &Plain) -> Result<Option<Arc<BaseGrid>>, Error> {
        Ok(ctx.memory_grids.get(name).cloned())
    }
}

pub(super) struct FileGridStore;
impl GridStore for FileGridStore {
    fn get_grid(&self, name: &str, ctx: &Plain) -> Result<Option<Arc<BaseGrid>>, Error> {
        let mut grids = GLOBALLY_ALLOCATED_GRIDS
            .get_or_init(init_grids)
            .lock()
            .unwrap();

        if let Some(grid) = grids.get(name) {
            return Ok(Some(grid));
        }

        let n = PathBuf::from(name);
        let ext = n
            .extension()
            .unwrap_or_default()
            .to_str()
            .unwrap_or_default();

        for path in &ctx.paths {
            let mut gridpath = path.clone();
            gridpath.push(name);
            let mut grid = std::fs::read(gridpath);

            if grid.is_err() {
                gridpath = path.clone();
                gridpath.push(ext);
                gridpath.push(name);
                grid = std::fs::read(gridpath);
            }
            let Ok(grid) = grid else {
                continue;
            };

            load_grid_into_collection(&mut grids, name, &grid)?;
            return Ok(grids.get(name));
        }

        Ok(None)
    }
}

pub(super) struct UnigridStore;
impl GridStore for UnigridStore {
    fn get_grid(&self, name: &str, ctx: &Plain) -> Result<Option<Arc<BaseGrid>>, Error> {
        for unigrid in &ctx.unigrid_elements {
            if let Some(grid) = unigrid.get(name) {
                return Ok(Some(grid.clone()));
            }
        }
        Ok(None)
    }
}

fn load_grid_into_collection(
    grids: &mut GridCollection,
    name: &str,
    bytes: &[u8],
) -> Result<(), Error> {
    let path = PathBuf::from(name);
    let ext = path.extension().and_then(|ext| ext.to_str());
    let grid = BaseGrid::read_bytes(name, bytes, ext)?;
    grids.0.insert(name.to_string(), Arc::new(grid));
    Ok(())
}
