
import { Route, Routes } from 'react-router-dom';
import MoleculeInfo from '../pages/Molecule';


  export const AppRoutes = ({ pages }) => (
    <Routes>
      <Route path='/' element={pages['Home']}></Route>
      {Object.keys(pages).map((page, index) => (
        <Route key={index} path={`/${page.toLowerCase()}`} element={pages[page]} />
      ))}
    </Routes>
  );

  export const MolRoutes = () => {
    return (
      <Routes>
        <Route path="molecule/:molid" element={<MoleculeInfo />} />
      </Routes>
    )
  }