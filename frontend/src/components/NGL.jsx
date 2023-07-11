import React, { useEffect, useRef, useState, useContext } from 'react';
import { Stage } from 'ngl';

const StageContext = React.createContext();

function NGLStage({ width, height, children }) {
  const stageRef = useRef(null);
  const [stage, setStage] = useState(null);

  useEffect(() => {
    if (stageRef.current) {
      const newStage = new Stage(stageRef.current, { backgroundColor: 'white' });
      setStage(newStage);

      return () => {
        newStage.dispose();
      };
    }
  }, []);

  return (
    <StageContext.Provider value={stage}>
      <div ref={stageRef} style={{ width, height }} />
      {children}
    </StageContext.Provider>
  );
}

function Component({ path, representationType = 'ball+stick' }) {
  const stage = useContext(StageContext);
  const componentRef = useRef(null);

  useEffect(() => {
    if (!stage) return;

    stage.loadFile(path).then(component => {
      componentRef.current = component;
      component.addRepresentation(representationType);
      stage.autoView();
    });

    return () => {
      if (componentRef.current) {
        stage.removeComponent(componentRef.current);
      }
    };
  }, [stage, path, representationType]);

  return null;
}

export {NGLStage, Component}